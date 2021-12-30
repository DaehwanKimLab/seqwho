#!/usr/bin/env python3
# --------------------------------------------------------------------------- #
# Built in Python 3.7.4                                                       #
# Copyright 2019, Christopher Bennett                                         #
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
# See <http://www.gnu.org/licenses/> for a copy of the GNU                    #
# General Public License                                                      #
# --------------------------------------------------------------------------- #

import os
import io
import re
import sys
import gzip
import glob
import copy
import json
import joblib
import random
import tarfile
import datetime
import subprocess
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as mplt
import multiprocessing as mp
from io import BytesIO
from argparse import ArgumentParser
from itertools import islice, product
from sklearn.ensemble import RandomForestClassifier

# Set Debug to true if searching for errors
debug = False
stat_nt2num = {'A' : 0, 
               'C' : 1, 
               'G' : 2, 
               'T' : 3, 
               'N' : 4}
# --------------------------------------------------------------------------- #
# Functions for detecting tags or overrepresented sequences in a file         #
# --------------------------------------------------------------------------- #
def build_posidx(k, nts = 'ACGT'):
    perm   = product(nts, repeat=k)
    nt2num = {}
    for i in range(len(nts)):
        nt2num[nts[i]] = i

    def pat2num(nt):
        if len(nt) == 0:
            return 0
        
        val = len(nts) * pat2num(nt[:-1]) + nt2num[nt[-1]]
        return val

    posidx = {}
    for muts in perm:
        mut = ''.join(muts)
        num = pat2num(mut)
        posidx[mut] = num

    return posidx

def tag_search(hashls, hashix, string):
    if not hashls or not hashix:
        hashls = [0 for _ in range(5 ** 6)]
        hashix = build_posidx(6, "ACGTN")

    for i in range(len(string)-6):
        sixmer = string[i:i+6]

        hashls[hashix[sixmer]] += 1

    return hashls, hashix

# --------------------------------------------------------------------------- #
# Set of functions for basic read handling and vector population              #
# --------------------------------------------------------------------------- #
def populate_vectors(all_files,         # List of files to type
                     index,             # Pass index into function
                     numlines = 100000, # Number of lines to read for index
                     building = False): # Read 40x more reads if not building  
    cnt   = 1
    stats = {}
    for fn in all_files:
        stats[fn] = {}
        # Stats
        stats[fn]["Per Base Seq"]   = seqpbase = []
        stats[fn]["Read lengths"]   = read_len = []
        stats[fn]["Quality Scores"] = quals    = []
        stats[fn]["Quality Dist"]   = qualcnt  = []
        polAhit  = 0
        readwNs  = 0
        readtot  = 0
        readpass = 0        
        hashix   = {}
        hashls   = []
        if not index.check_file(fn):
            print(">>>> skipping file %s" % (fn), file=sys.stderr)
            continue 

        print("> Running file %d: %s" % (cnt, fn), file=sys.stderr)
        cnt += 1

        snag     = numlines * 40 if not building else numlines
        isread   = False
        isqual   = False
        isfastq  = False
        everhit  = False
        filesize = os.path.getsize(fn)
        index.set_working_file(fn)

        ific = open(fn, "rb")
        ifi  = gzip.GzipFile(fileobj=ific)
        with io.TextIOWrapper(ifi, encoding='utf-8') as ifo:
            chunk   = list(islice(ifo, snag))
            linecnt = 0
            for line in chunk:
                line    = line.strip().strip('\n')
                linelen = len(line)
                if linelen == 0:
                    continue

                linecnt += 1

                if isqual and not building:
                    if len(quals) < linelen: # expand the quals if needed
                        quals += [0 for _ in range(linelen - len(quals))]
                    for i in range(linelen):
                        qual = ord(line[i])-33
                        if len(qualcnt) < qual:
                            qualcnt += [0 for _ in range(qual - len(qualcnt))]
                        elif len(qualcnt) == 0: # fixes issue when first qual is 0
                            qualcnt += [0]
                        quals[i] += qual
                        qualcnt[qual-1] += 1

                    isqual  = False
                    isfastq = True
                    continue

                if isread:
                    if re.search('[^ACGTN]', line) or linelen < 10:
                        readpass += 1
                        continue
                    
                    # Update stats
                    if not building:
                        readtot += 1   
                        if 'A' * linelen == line:
                            polAhit += 1 
                        if "N" in line:
                            readwNs += 1 
                        if len(seqpbase) < linelen:
                            seqpbase += [[
                                          0, # A
                                          0, # C
                                          0, # G
                                          0, # T
                                          0, # N
                                          0  # Total
                                         ] for _ in range(linelen - len(seqpbase))]
                        
                        for i in range(linelen):
                            seqpbase[i][stat_nt2num[line[i]]] += 1
                            seqpbase[i][5] += 1
                        
                        if len(read_len) < linelen:
                            read_len += [0 for _ in range(linelen - len(read_len))]
                            
                        read_len[linelen - 1] += 1
                        hashls, hashix = tag_search(hashls, hashix, line[:10])

                    # Update index
                    if linecnt < numlines:
                        if len(line) > 50:
                            line = line[10:-10]
                        else:
                            line = line[10:]
                        index.add_read_to_index(line)
                    isread  = False
                    everhit = True
                    continue

                if line[0] in [">", "@"]:
                    isread = True
                if line.startswith("+"):
                    isqual = True

            if not building:
                file_bytes = ific.tell() if len(chunk) == snag else filesize
                nread = int(len(chunk)/4) if isfastq else int(len(chunk)/2)
                stats[fn]['Estimated Read Number'] = int(nread 
                                                          * filesize 
                                                          / file_bytes)
                stats[fn]['File format'] = "fastq" if isfastq else "fasta"

        index.close_file(fn)
        ific.close()
        ifi.close()

        if not everhit:
            stats[fn] = {"Error" : "Possible Colorspace read"}
            continue

        if not building:
            tag_report = []
            tagsum     = sum(hashls)
            for tag, i in hashix.items():
                freq = round(hashls[i] / tagsum, 4)
                if freq > 0.1:
                    tag_report.append((tag, freq))
            stats[fn]["Biased 5' end 6-mers"] = tag_report

            sum_lens  = 0
            sum_quals = 0
            tot_bases = 0
            gc_bases  = 0
            for i in range(len(seqpbase)):
                tot        = seqpbase[i][5]
                sum_lens  += (i + 1) * read_len[i]
                sum_quals += quals[i]
                tot_bases += tot
                quals[i]   = round(quals[i]/tot, 4)
                for j in range(5):
                    if j in (1, 2): # Looking for C and G
                        gc_bases += seqpbase[i][j]
                    seqpbase[i][j] = round(seqpbase[i][j] / tot, 4)

            stats[fn]["Mean Quality"]   = round(sum_quals / tot_bases, 3)
            stats[fn]["Mean Read Len"]  = round(sum_lens  / readtot,   3)
            stats[fn]["Perc PolyA"]     = round(polAhit   / readtot,   5) * 100
            stats[fn]["Perc GC Cont"]   = round(gc_bases  / tot_bases, 5) * 100
            stats[fn]["Perc Reads w N"] = round(readwNs   / readtot,   5) * 100
            stats[fn]["Reads Passed"]   = readtot
            stats[fn]["Reads Omitted"]  = readpass
            stats[fn]["Perc Passed"]    = round(readtot / (readpass + readtot), 5) * 100
        
    return index, stats

# --------------------------------------------------------------------------- #
# Index Vector Classes                                                        #
#                                                                             #
# Will build a base array to load 1-5mer frequencies into for each file       #
# --------------------------------------------------------------------------- #
class BaseVector(object):
    # Internal required variables variables:
    # k      = max kmer length
    # nts    = Nucleotide Alphabet
    # nt2num = Conversion from Char to Int
    # arr    = Vector or Array with positions linked to String
    # nt2pos = String or Char map to Integer positon in arr
    # bounds = ranges in vector corresponding to data types

    #--------------------------------------- #
    # Public Functions for class             #
    #--------------------------------------- #
    def initialize_vector(self, k_size, nt_list):
        self.k   = k_size
        self.nts = nt_list

        self._build_nt2num(nt_list)
        self._build_arr()
        self._build_posidx()

        self.bounds = [[0, len(self.arr)]]

    def add_to_ix(self, kmer):
        if kmer not in self.nt2pos:
            return

        self.arr[self.nt2pos[kmer]] += 1

    def set_ix(self, normalize = True):
        if normalize:
            return self._normalize_ix()

        return self.arr[:]

    #--------------------------------------- #
    # Hidden Functions for class             #
    #--------------------------------------- #
    def _normalize_ix(self):
        first = 0
        last  = 0
        for i in range(self.k):
            # kmer size X length of alphabet + first index = last value of set
            last = (len(self.nts) ** (i + 1)) + first 
            tot  = sum(self.arr[first:last])
            if tot == 0:
                return [0]

            while first < last:
                self.arr[first] = self.arr[first]/tot
                first += 1

        return self.arr[:]

    # Convert a k-mer string to number and find the position 
    # it would have in a concatinated vector of smaller kmers
    def _pat2num(self, 
                 nt, 
                 tnt = None, 
                 klen = None):
        if len(nt) == 0:
            return 0
        
        if tnt is None:             # This is only run at the upper 
            tnt  = len(self.nts)     # recursion of the function
            klen = len(nt) - 1
        
        val = tnt * self._pat2num(nt[:-1], tnt) + self.nt2num[nt[-1]]

        if klen is None:            # This should only be none in 
            return val              # recursive loops of pat2num
        
        while klen > 0:             # Find Position in concatinated 
            val  += tnt ** klen      # vector of smaller kmers 
            klen -= 1

        return val   

    def _build_nt2num(self, nts = 'ACGT'):
        nt2num = {}
        for i in range(len(nts)):
            nt2num[nts[i]] = i
        
        self.nt2num = nt2num

    def _build_arr(self):
        k   = self.k
        l   = len(self.nts)
        arr = []
        while k > 0:
            arr += [0 for x in range(l ** k)]
            k   -= 1

        self.arr = arr

    def _build_posidx(self):
        k      = self.k
        posidx = {}        
        while k > 0:
            perm = product(self.nts, repeat=k)
            for muts in perm:
                mut = ''.join(muts)
                num = self._pat2num(mut)
                posidx[mut] = num
            k -= 1

        self.nt2pos = posidx 

# --------------------------------------------------------------------------- #
# Used to build repeat index                                                  #
#                                                                             #
# Only one class is used and will store information for each file             #
# --------------------------------------------------------------------------- #
class RepVector(object):
    # Internal required variables variables:
    index     = {}
    buildmode = False

    #--------------------------------------- #
    # Public Functions for class             #
    #--------------------------------------- #

    # Build index and set buildmode to true to allow purging of the data set.
    # This also requires a set of known labels
    def initialize_build_vector(self, repeats, labelfn):
        self.flabels, \
          self.labkey, \
          self.labelmetadata, \
          self.spe_cnt, \
          self.typ_cnt \
            = self._get_labels(labelfn)
        
        self.repeats    = repeats
        self.buildmode  = True

    # Initialize vector using a prebuilt index
    def initialize_vector_from_index(self, index, fnlist):
        self.flabels        = set(fnlist)
        self.labkey         = index["label key"]
        self.labelmetadata  = index["label metadata"]

        self.spe_cnt        = index["label metadata"]["nspecies"]
        self.typ_cnt        = index["label metadata"]["ntypes"]
        self.k              = len(index["kmer sizes"])
        self.repsize        = index["repeat size"]
        self.repeats        = set()

        for repeat in index["repeat list"]:
            self.index[repeat] = {}

    def get_fnames(self):
        names = []
        for fn in self.flabels:
            names.append(fn)

        return(names)

    def check_file(self, fn):
        if fn in self.flabels:
            return True

        return False

    def add_to_reps(self, 
                    kmer, 
                    fn, 
                    cnt = 0):
        if kmer not in self.index:
            if kmer not in self.repeats:
                return
            
            self.repeats.remove(kmer)
            self._initialize_entry(kmer, [0 for x in range(self.spe_cnt)])

        if fn not in self.index[kmer]:
            self._initialize_entry(kmer, fn)

        self.index[kmer][fn] += 1

        if self.buildmode:
            spec, seqtype = self.flabels[fn]
            ix = self.labkey[spec]
            self.index[kmer]['total'][ix] += 1

            if cnt % 100 == 0 and cnt != 0:
                self._purge_kmers(cnt)

    #--------------------------------------- #
    # Hidden Functions for class             #
    #--------------------------------------- #
    def _purge_kmers(self, cnt):
        rmmer_ = set()
        for k, counts in self.index.items():
            if sum(counts['total']) < cnt:
                rmmer_.add(k)

        while rmmer_:
            rm = rmmer_.pop()
            del self.index[rm]

    def _initialize_entry(self, kmer, entry):
        if isinstance(entry, list):
            self.index[kmer] = {'total' : copy.deepcopy(entry)}
        else:
            self.index[kmer][entry] = 0

    # Extract labels from a list of files and labels in format: 
    # $file_ID,$species,$seqtype
    @staticmethod
    def _get_labels(fname):
        if not os.path.exists(fname):
            print("Cannot find file %s" % fname)
            exit(0)

        labs = open(fname, "r").read()
        labs = labs.strip().split("\n")

        ## CD - debug
        # labs = random.choices(labs, k=200)

        key_collector = {"species" : [], "type" : []}
        flabels       = {}
        spe_cnt       = 0
        typ_cnt       = 0
        for lab in labs:
            f, spec, seqtype = lab.split(",")
            if spec not in key_collector["species"]:
                key_collector["species"].append(spec)
                spe_cnt += 1 
                
            if seqtype not in key_collector["type"]:
                key_collector["type"].append(seqtype)
                typ_cnt += 1
                
            flabels[f + ".fastq.gz"] = (spec, seqtype)

        labkey = {}
        ixcnt  = 0
        for category, sets in key_collector.items():
            for s in sets:
                if s not in labkey:
                    labkey[s] = ixcnt
                    ixcnt += 1

        key_collector['nspecies'] = spe_cnt 
        key_collector['ntypes']   = typ_cnt

        return flabels, labkey, key_collector, spe_cnt, typ_cnt

    # Normalize repeat index by scaling the data
    @staticmethod
    def _normalize_ix(arr):
        max_int = max(arr)
        min_int = min(arr)
        if max_int != 0 or (max_int-min_int) != 0:
            for i in range(len(arr)):
                arr[i] = (arr[i] - min_int)/(max_int - min_int) # scale data
        
        return arr

class IxVector(RepVector):
    #--------------------------------------- #
    # Public Functions for class             #
    #--------------------------------------- #
    def init_new(self, k, ks, repeats, labelfn, useN = False):
        self.alph       = 'ACGTN' if useN else 'ACGT'
        self.k          = k
        self.repsize    = ks
        self.sizes      = self._get_sizes(k,ks)

        self.basevect   = BaseVector()
        self.basevect.initialize_vector(self.k, self.alph)

        repeatset       = self._kmerize_repeats(repeats, ks)
        self.initialize_build_vector(repeatset, labelfn)

        self.fns        = {}
        self.k_ix       = []
        self.working_ix = None
        self.working_fn = None

    def load_index_params(self, index, fnlist, useN = False):
        self.alph       = 'ACGTN' if useN else 'ACGT'
        self.fns        = {} 
        self.k_ix       = []
        self.working_ix = None
        self.working_fn = None

        self.initialize_vector_from_index(index, fnlist)
        self.basevect = BaseVector()

        self.basevect.initialize_vector(self.k, self.alph)
        self.sizes = self._get_sizes(self.k, self.repsize)

    def set_working_file(self, fn):
        if fn in self.fns:
            self.working_ix = self.fns[fn]
            self.working_fn = fn
            return

        self.working_ix = len(self.k_ix)
        self.working_fn = fn

        self.fns[fn] = self.working_ix
        self.k_ix.append(copy.deepcopy(self.basevect))

    def close_file(self, fn):
        ix = self.fns[fn]
        self.k_ix[ix] = self.k_ix[ix].set_ix()

    def add_read_to_index(self, read):
        if self.working_ix is None:
            print("Working file not initialized with set_working_file", 
                  file=sys.stderr)
            exit(0)
            
        ix = self.working_ix 
        for i in range(len(read)):
            for j in self.sizes:
                if i + j > len(read):
                    continue

                kmer = read[i:i+j]

                if len(kmer) == self.repsize:
                    self.add_to_reps(kmer, self.working_fn, ix)
                elif len(kmer) <= self.k:
                    self.k_ix[ix].add_to_ix(kmer)
                else:
                    print("Impossible state entered", 
                          file=sys.stderr)
                    print("Error in size of kmers", 
                          file=sys.stderr)
                    exit(1)

    def set_ix(self, norm = True):
        matrix = [[0.0 for x in range(len(self.index))] \
                    for y in range(len(self.fns))]
        repix = sorted([key for key in self.index])
        
        for i in range(len(repix)):
            rep = repix[i]
            for fn, count in self.index[rep].items():
                if fn == 'total':
                    continue

                matrix[self.fns[fn]][i] += count

        self.index.clear()

        labs  = [] 
        vects = []
        for fn, ix in self.fns.items():
            if self.buildmode:
                labs.append([self.labkey[x] for x in self.flabels[fn]])

            if norm:
                mx = self._normalize_ix(matrix[ix])
            else:
                mx = matrix[ix]

            kmx = self.k_ix[ix]

            vects.append(kmx + mx)

        if self.buildmode:
            basic_ix = {
                "kmer sizes"     : self.sizes[:-1],
                "repeat size"    : self.repsize, 
                "label key"      : self.labkey,  
                "label metadata" : self.labelmetadata,
                "repeat list"    : repix,
                "file ix"        : self.fns
            }

            return (np.array(labs), 
                    np.array([np.array(vec) for vec in vects]), 
                    basic_ix)
        
        else:
            return (self.fns, 
                    np.array([np.array(vec) for vec in vects]))

    def purge_kmers(self):
        if not self.buildmode:
            print("Cannot purge kmers when not in build mode", 
                  file=sys.stderr)
            return

        nspecies = self.spe_cnt

        def variance(vector):
            n       = len(vector)
            mean    = sum(vector)/n
            sigmasq = 0
            for v in vector:
                sigmasq += (v-mean)**2

            sigmasq = sigmasq/n
            return sigmasq

        perfect    = [0.0 for x in range(nspecies)]
        perfect[0] = 1.0
        maxvar = variance(perfect)
        cutoff = maxvar / 2    
        kept_k = {}
        for kmer, counts in self.index.items():
            total = sum(counts['total'])   
            if total == 0 or total < 50:
                continue

            if kmer == 'A' * len(kmer):
                continue

            freqs = [(x/total) for x in counts['total']]

            v = variance(freqs)
            if v > cutoff:
                kept_k[kmer] = copy.deepcopy(counts)
                kept_k[kmer]['total'] = freqs + [v]

        self.index.clear()
        self.index = kept_k

    def check_state(self):
        print("In build mode?")
        return self.buildmode

    #--------------------------------------- #
    # Hidden Functions for class             #
    #--------------------------------------- #
    @staticmethod
    def _get_sizes(k, repsize):
        sizes = []
        for i in range(k):
            sizes.append(i+1)
        
        sizes.append(repsize)
        return(tuple(sizes))

    @staticmethod
    def _kmerize_repeats(repeats, k_size):
        # Read and count kmers in repeat index
        kmers = set()

        if debug: # CB-debug
            ofo = open("kmer_array_prefilter.empty", "w")

        while repeats:
            read = repeats.pop()
            for k in range(len(read)-k_size+1):
                kmer = read[k:k+k_size]
                if kmer not in kmers:
                    kmers.add(kmer)
                    if debug: ofo.write(kmer + "\n")

        if debug: ofo.close()
        return kmers

    @staticmethod
    def _minimize(repeats, size):
        hashset = set()

        while repeats:
            rep        = repeats.pop()
            rep_size   = len(rep)
            rep_counts = {}
            for i in range(rep_size-size+1):
                kmer = rep[i:i+size]
                if kmer not in rep_counts:
                    rep_counts[kmer] = 0

                rep_counts[kmer] += 1

            max_int  = -1
            max_kmer = ''
            for k, c in rep_counts.items():
                if c > max_int:
                    max_int = c
                    max_kmer = k

            hashset.add(max_kmer)

        return hashset    

# --------------------------------------------------------------------------- #
# Class for running the tests and extracting statistics for each file         #
#                                                                             #
# Will wrap the other classes for building and set internal state to false    #
# --------------------------------------------------------------------------- #
class SeqWho_Index(object):
    def __init__(self, fnindex):
        # keys is the json index that contains all metadata and repeat identity
        # random forests contains all random forest models
        tars = tarfile.open(fnindex, "r:gz")
        keys = tars.extractfile("index.json")
        self.stats   = {}
        self.keys    = json.load(keys)
        self.convkey = self._dict_swap(self.keys["label key"])
        
        nrf = self.keys['model file']['max models']
        self.random_forests      = [None for x in range(nrf)]
        self.full_random_forests = [None, None]

        for name, model_info in self.keys["model file"].items():
            if name == "max models":
                continue

            ix, modfile = model_info
            f = tars.extractfile(modfile)            
            if "species" in name and "full_set" in modfile:
                self.full_random_forests[0] = joblib.load(f)
            elif "type" in name and "full_set" in modfile:
                self.full_random_forests[1] = joblib.load(f)
            else:
                self.random_forests[ix] = joblib.load(f)

        tars.close()

        self.vectors = None

    #--------------------------------------- #
    # Public Functions for class             #
    #--------------------------------------- #
    def init_vectors(self, fnames):
        vects = IxVector()
        vects.load_index_params(self.keys, fnames)

        self.vectors = vects

    def run_vector_population(self):
        numlines = 100000
        if "NumberLine" in self.keys:
            numlines = self.keys["NumberLine"]
        # try:
        self.vectors, self.stats = populate_vectors(self.vectors.get_fnames(),
                                                    self.vectors,
                                                    numlines)
        #     return True
        # except Exception as exc:
        #     print("Unexpected Error:", exc,
        #          file=sys.stderr)
        #     pass

        return True

    def run_vector_test(self, outfn):
        fns, mtx = self.vectors.set_ix()
        assert len(fns) == len(mtx)

        callresults = []
        for fn, i in fns.items():
            if "Error" in self.stats[fn]:
                result = self._get_file_table(fn)
                callresults.append(result)
                continue

            newadd     = []
            confidence = "high"
            for j in range(len(self.random_forests)):
                model = self.random_forests[j]

                try:
                    res = model.predict(mtx[i].reshape(1, -1))
                    if res[0] == 1:
                        newadd.append(self.convkey[str(j)])
                except:
                    newadd.append("list_error")
                    break

            if len(newadd) != 2:
                wgt = [0 for x in newadd]
                for j in range(len(self.full_random_forests)):
                    model = self.full_random_forests[j]

                    try:
                        res = model.predict(mtx[i].reshape(1, -1))
                        res = int(res[0])
                        if j == 1:
                            res += self.keys["label metadata"]["nspecies"]
                        
                        callpart = self.convkey[str(res)]
                        if callpart not in newadd:
                            newadd.append(callpart)
                            wgt.append(0)
                            confidence = "low:"
                        else:
                            wgt[newadd.index(callpart)] += 1

                    except:
                        newadd.append("list_error")
                        confidence = "NA"
                        break

                if confidence != "NA":
                    confidence = confidence + ":".join([str(x) for x in wgt])

            self._add_call_to_stats(fn, newadd, confidence)
            result = self._get_file_table(fn)
            callresults.append(result)

        self._plot_stats(outfn)
        self._save_json(outfn)
        self._save_table(callresults, outfn)
        #return callresults

    #--------------------------------------- #
    # Hidden Functions for class             #
    #--------------------------------------- #
    @staticmethod
    def _dict_swap(d):
        new_d = {}
        for key, val in d.items():
            new_d[str(val)] = key
        
        return new_d

    def _save_json(self, out):
        with open("%s.json" % out, "w") as ofo:
            json.dump(self.stats, ofo, indent=1)     

    def _save_table(self, table, out):
         with open("%s.tsv" % out, "w") as ofo:
             for line in table:
                 ofo.write(line + "\n")   

    def _plot_stats(self, out):
        if not self.stats:
            print("Error: no data to plot",
                  file=sys.stderr)
            exit(1)

        outpath = out + "_plots"
        if not os.path.exists(outpath):
            os.mkdir(outpath)

        for fn, stat in self.stats.items():
            if "Error" in stat:
                continue

            if "/" in fn:
                fn = fn.split("/")[-1]

            sb.set_style('whitegrid')
            sb.set(font_scale=3)
            fig, (ax1, ax2, ax3, ax4) = mplt.subplots(4,1,figsize=(25,30))
            qualdist = {"Count"         : np.array(stat["Quality Dist"]),
                        "Quality Score" : np.array(
                            [i for i in range(len(stat["Quality Dist"]))]
                        )}
            qualdist = pd.DataFrame(qualdist)
            sb.lineplot(x = "Quality Score", 
                        y = "Count", 
                        data = qualdist, 
                        ax=ax1,
                        linewidth=4,
                        color = "blue")

            qualscore = {"Quality Score" : np.array(stat["Quality Scores"]),
                         "Pos"           : np.array(
                             [i for i in range(len(stat["Quality Scores"]))]
                         )}
            qualscore = pd.DataFrame(qualscore)
            sb.lineplot(x = "Pos", 
                        y = "Quality Score", 
                        data = qualscore, 
                        ax=ax2,
                        linewidth=4,
                        color = "blue")
            ax2.set(ylim=(0, 40))

            readdist = {"Count"  : np.array(stat["Read lengths"]),
                        "Length" : np.array(
                            [i for i in range(len(stat["Read lengths"]))]
                        )}
            readdist = pd.DataFrame(readdist)
            sb.lineplot(x = "Length", 
                        y = "Count", 
                        data = readdist, 
                        ax=ax3,
                        linewidth=4,
                        color = "blue")

            seqdist = pd.DataFrame(stat["Per Base Seq"]) 
            seqdist.columns = ['A', 'C', 'G', 'T', 'N', 'Total']
            seqdist = seqdist.drop('Total', axis = 1).transpose()
            seqdist.columns = [str(i) for i in range(len(stat["Per Base Seq"]))]
            sb.heatmap(seqdist, 
                       cmap="YlGnBu", 
                       vmin=0, 
                       vmax=1, 
                       ax=ax4, 
                       cbar_kws = {'label' : "Prevelence"})
            mplt.tight_layout()
            mplt.savefig(outpath + "/" + fn + "_statistics.png")
            mplt.clf() # Clear Figure
            mplt.close('all')

            ## TODO add method of combining plots into same axis

    def _get_file_table(self, fn):
        if fn not in self.stats:
            print("Error, no file %s found in stats" % fn,
                  file=sys.stderr)
            exit(1)
        
        # Convert everything at level to string
        def stringify(name, level):
            string = '%s: ' % name
            try:
                coll = []
                for key, val in level.items():
                    if isinstance(val, list):
                        val = "-".join(val)
                    coll.append("%s %s: %s" % (name, key, val))
                return '\t'.join(coll)
            except Exception as exc:
                # print(exc)
                pass
            if isinstance(level, list):
                try:
                    string += ";".join([
                                    "-".join([str(l) for l in lev])
                                for lev in level])
                    return string
                except Exception as exc:
                    # print(exc)
                    pass
                try:
                    string += ";".join([str(lev) for lev in level])
                    return string
                except Exception as exc:
                    # print(exc)
                    pass
            string += "%s" % str(level)
            return string

        stat = self.stats[fn]
        report = []
        for item, value in stat.items():
            report.append(stringify(item, value))

        report = "\t".join(report)
        return fn + "\t" + report

    def _add_call_to_stats(self, fn, call, confidence):
        if "Call" in self.stats[fn]:
            print("Error call already associated with a file",
                  file=sys.stderr)
            exit(1)

        self.stats[fn]['Call'] = {"Confidence" : confidence,
                                  "Species"    : [],
                                  "Seq Type"   : []}
        for c in call:
            if c == "list_error":
                self.stats[fn]["Call"] = "No Call"
            elif c in self.keys["label metadata"]["species"]:
                self.stats[fn]['Call']["Species"].append(c)
            else:
                assert c in self.keys["label metadata"]["type"], "%s call not found" % c
                self.stats[fn]['Call']["Seq Type"].append(c)


