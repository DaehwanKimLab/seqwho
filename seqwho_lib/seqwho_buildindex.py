#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
# Built in Python 3.7.4                                                       #
# Copyright 2019, Christopher Bennett                                         #
#                                                                             #
# This script build the SeqWho index needed for typing sequencing files       #
#                                                                             #
# Seq-Who is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# Seq-Who is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
# GNU General Public License for more details.                                #
#                                                                             #
# See <http://www.gnu.org/licenses/> for a copy of the GNU General            #
# Public License                                                              #
# --------------------------------------------------------------------------- #

import sys
import os
import subprocess
import re
import gzip
import glob
import copy
import random
import json
import tarfile
import datetime
import joblib
import multiprocessing as mp
import pandas as pd
import numpy as np
from itertools import islice
from itertools import product
from argparse import ArgumentParser
from io import BytesIO
from sklearn.ensemble import RandomForestClassifier
import seqwho_modules as sqwho

# --------------------------------------------------------------------------- #
# Base common functions                                                       #
# --------------------------------------------------------------------------- #
def save_matrix(labels, 
                matrix, 
                keys, 
                fn):
    fnmtx = fn + "_mtx.csv"
    fnlab = fn + "_lab.csv"
    fnkey = fn + "_key.json"

    with open(fnkey, "w") as ofo:
        json.dump(keys, ofo, indent=1)

    mx_ofo = open(fnmtx, "w")
    la_ofo = open(fnlab, "w")
    for i in range(len(labels)):
        mx = matrix[i]
        label = labels[i]

        mtxline = "\t".join([str(x) for x in mx])
        labline = "\t".join([str(x) for x in label])

        mx_ofo.write(mtxline + '\n')
        la_ofo.write(labline + '\n')

    mx_ofo.close()
    la_ofo.close()

    return([fnkey, fnlab, fnmtx])

# --------------------------------------------------------------------------- #
# Main functions for extracting repeats                                       #
# --------------------------------------------------------------------------- #
def extract_repeats(h2_ridx):
    file_out = h2_ridx.replace(".seed", ".txt")
    repeat_extract = "sed 's/\x1b\[[0-9;]*m//g' %s | \
                      awk '{OFS=\"\t\"}{if ($5 == \"+\") \
                          {print \"chr\"$10,$11,$12,\".\",$5} \
                        else {print \"chr\"$11,$10,$12,\".\",$5}}' | \
                      awk '{OFS=\"\t\"} \
                          sub(/\:/,\"\t\"){$1=$1}1 \
                          sub(/\:/,\"\t\"){$2=$2}2' | \
                      cut -f5 > %s" % (h2_ridx, file_out)
    os.system(repeat_extract)
    if not os.path.exists(file_out):
        print("Failed to Build")
        exit(1)

def load_repeats(h2_ridx):
    repeats = set()
    for repeatsfn in h2_ridx:
        rep = open(repeatsfn, "r").read()
        rep = rep.split("\n")
        repeats.update(rep)

    return repeats

# --------------------------------------------------------------------------- #
# Main Functions for building Index                                           #
# --------------------------------------------------------------------------- #
def build_input(h2_ridx, 
                labfile, 
                mask, 
                k_size, 
                repk,
                numlines,
                outfn):
    # Build Repeats and master index
    print("Building")
    h2_ridx = h2_ridx.split(",")
    if not isinstance(h2_ridx, list):
        h2_ridx = [h2_ridx]

    for repeatsfn in h2_ridx:
        if repeatsfn.endswith(".seed"):
            if not os.path.exists(repeatsfn.replace('seed', 'txt')):
                extract_repeats(repeatsfn)
            else:
                repeatsfn = repeatsfn.replace('seed', 'txt')

    rep_sets = load_repeats(h2_ridx)

    index = sqwho.IxVector()
    index.init_new(k_size, 
                   repk, 
                   rep_sets, 
                   labfile)

    # Make list of file index links
    all_files = glob.glob("*fastq.gz")
    random.shuffle(all_files)

    index, stats = sqwho.populate_vectors(all_files, index, numlines, building=True)
    index.purge_kmers()

    labels, matrix, keys = index.set_ix()
    keys["index name"] = outfn
    keys["NumberLine"] = numlines

    return(save_matrix(labels, 
                       matrix, 
                       keys, 
                       "intermediate_" + outfn + "/kbloom_vectors"))

# --------------------------------------------------------------------------- #
# Build the Random Forest models for each selection                           #
# --------------------------------------------------------------------------- #
def build_model(repfn,
                labelfn,
                mask,
                ksize,
                repsize,
                numlines,
                outfn):
    keyfn, labelfn, matrixfn = build_input(repfn,
                                           labelfn,
                                           mask,
                                           ksize,
                                           repsize,
                                           numlines,
                                           outfn)

    ## CB - DEBUG - If files are already built using default file names
    ## this can be used to test random forest build. Comment out the
    ## above code and remove the # on the three lines below
    # keyfn = "intermediate_SeqWho.ix/kbloom_vectors_key.json"
    # labelfn = "intermediate_SeqWho.ix/kbloom_vectors_lab.csv"
    # matrixfn = "intermediate_SeqWho.ix/kbloom_vectors_mtx.csv"
    if not os.path.exists(keyfn):
        print("No key file", file=sys.stderr)
        exit(1)
    elif not os.path.exists(labelfn):
        print("No label file", file=sys.stderr)
        exit(1)
    elif not os.path.exists(matrixfn):
        print("No matrix file", file=sys.stderr)
        exit(1)

    with open(keyfn, "r") as ifi:
        keyset = json.load(ifi)

    sets    = pd.read_csv(matrixfn, sep="\t", header=None, engine="python")
    alabels = pd.read_csv(labelfn, sep="\t", header=None)
    mask    = sets.notna().all(axis=1)
    sets    = sets[mask]
    alabels = alabels[mask]
    tar_out = tarfile.open(outfn, "w:gz")

    # Add models to tar file index as pickles
    keyset["model file"] = {}
    prediction_tester    = {}
    model_index = 0
    for lab, keys in keyset["label metadata"].items():
        if lab == "species":
            labels      = alabels.iloc[:,0]
            full_labset = labels
        elif lab == "type":
            labels = alabels.iloc[:,1]

            ## subtract nspecies so we have a full correct label set at min zero
            full_labset = labels - keyset["label metadata"]["nspecies"]                
        elif lab in ["nspecies", "ntypes"]:
            continue
        else:
            print("Invalid index! Try rebuilding",
                  file=sys.stderr)
            exit(1)

        ## Save full Random forest model for the categories 
        ## (a double check method)
        filename = "full_set_" + lab + "_rf_mod.pkl"
        clf      = RandomForestClassifier(random_state=0, 
                                          n_jobs=-1, 
                                          n_estimators=500)
        model    = clf.fit(sets, full_labset)

        keyset["model file"][lab] = [model_index, filename]
        bytes_container = BytesIO()
        joblib.dump(model, bytes_container)
        bytes_container.seek(0)   
        info       = tarfile.TarInfo(name=filename)
        info.mtime = datetime.datetime.timestamp(datetime.datetime.now())
        info.size  = bytes_container.getbuffer().nbytes
        tar_out.addfile(tarinfo=info,fileobj=bytes_container)

        ## Save individual random forest models for categories
        for name in keys:
            categ    = keyset["label key"][name]
            filename = name + "_rf_mod.pkl"

            keyset["model file"][name]         = [model_index, filename]
            keyset["model file"]["max models"] = model_index + 1
            
            bitlab = np.where(labels == categ, 1, 0)
            clf    = RandomForestClassifier(random_state=0, 
                                            n_jobs = -1, 
                                            n_estimators=500)
            model  = clf.fit(sets, bitlab)
            
            bytes_container = BytesIO()
            joblib.dump(model, bytes_container)
            bytes_container.seek(0)  
            info       = tarfile.TarInfo(name=filename)
            info.mtime = datetime.datetime.timestamp(datetime.datetime.now())
            info.size  = bytes_container.getbuffer().nbytes
            tar_out.addfile(tarinfo=info,fileobj=bytes_container)
            
            model_index += 1
            prediction_tester[filename] = model.predict(sets)            

    # Add Keys/Index to tar file index in json format
    bytes_container = BytesIO(json.dumps(keyset).encode())
    bytes_container.seek(0)
    info       = tarfile.TarInfo(name="index.json")
    info.mtime = datetime.datetime.timestamp(datetime.datetime.now())
    info.size  = bytes_container.getbuffer().nbytes
    tar_out.addfile(tarinfo=info,fileobj=bytes_container)
    
    tar_out.close()

    ## CB DEBUG- This block tests the fidelity of model saving and 
    ## ensures it is loaded faithfully
    # tars = tarfile.open(outfn, "r:gz")
    # keys = tars.extractfile("index.json")

    # keys = json.load(keys)
    # # print(keys["model file"])
    # for name, model_info in keys["model file"].items():
    #     if name == "max models":
    #         continue

    #     print(name, model_info)
    #     f = tars.extractfile(model_info[1])
    #     model = joblib.load(f)
    #     print(all(prediction_tester[model_info[1]] == model.predict(sets)))

    # tars.close()

def main():
    parser = ArgumentParser(
        description='Build SeqWho index for classifying sequence file types')
    parser.add_argument('-r', '--repeats',
                        dest = 'repfn',
                        type=str,
                        required=True,
                        help = 'Comma separated list of repeat indecies from \
                                HISAT2 for each species') 
    parser.add_argument('-l', '--labels',
                        dest = 'labelfn',
                        type=str,
                        required=True,
                        help = 'CSV file with: base_filename, species, \
                                sequence_file_type') 
    parser.add_argument('-m', '--mask',
                        dest = 'mask',
                        type=str,
                        help = 'Comma separated list of any file types \
                                you want to ommit')                                                 
    parser.add_argument('-k', '--ksize',
                        dest = 'ksize',
                        default=5,
                        help = 'Size of full kmer vector (default: 5)')
    parser.add_argument('-j', '--repeat-kmer',
                        dest = 'repsize',
                        default=31,
                        help = 'Size of repeat kmers to use (default: 31)') 
    parser.add_argument('-n', '--number_of_lines',
                        dest = 'numlines',
                        default=100000,
                        help = 'Number of lines to draw from the files. (Default: 100000')    
    parser.add_argument('-o', '--out',
                        dest="outfn",
                        default="SeqWho.ix",
                        type=str,
                        help="Final name of SeqWho index file")
    parser.add_argument('--rebuild',
                        dest="rebuild",
                        action="store_true",
                        help="If SeqWho index is detected build will stop.")
    parser.add_argument('-v', '--verbose',
                        dest = 'verbose',
                        action = "store_true",
                        help = 'Verbose Output') 

    args = parser.parse_args()

    outfile = "intermediate_" + args.outfn
    if os.path.exists(args.outfn) or os.path.exists(outfile):
        if args.rebuild:
            try:
                os.remove(args.outfn)
            except:
                pass

            try:
                assert outfile != "/"
                os.system("rm -rf %s" % outfile)
            except:
                pass

        else:
            print("SeqWho Index files detected. Exiting Build",
                  file=sys.stderr)
            exit(0)

    os.mkdir(outfile)

    # ---------------------------------------------------- #
    # Basic setup to build index and index building        #
    # ---------------------------------------------------- #
    if args.mask:
        mask = set(args.mask.split(","))
    else:
        mask = []
   
    build_model(args.repfn,
                args.labelfn,
                mask,
                args.ksize,
                args.repsize,
                args.numlines,
                args.outfn)

    return 0

if __name__ == '__main__':
    sys.exit(main())
