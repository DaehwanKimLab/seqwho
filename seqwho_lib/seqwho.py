#!/usr/bin/env python
# --------------------------------------------------------------------------- #
# Built in Python 3.7.4                                                       #
# Copyright 2019, Christopher Bennett                                         #
#                                                                             #
# The purpose of this script is to perform aligment on sequences from unknown #
# files for summary statistics and to attempt typing on them                  #                                                          
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

import sys
import os
import subprocess
import re
import random
import gzip
import glob
import multiprocessing
import argparse
import seqwho_modules as sqwho
from itertools import islice, product

cigar2num = {'S':0,'D':1,'I':2,'N':3}
# def alignment_stats(fname, 
#                     species, 
#                     aligner,
#                     paired):
def alignReads(genome, read_chunk, longest):
    cigar_re    = re.compile('\d+\w')
    metacigar   = [[0 for _ in 'SDIN'] for _ in range(longest)]
    cigarcounts = {}

    def extend_meta(ix, code, count):
        metacigar[ix][cigar2num[code]] += count

    if not glob.glob(genome + "*"):
        print("%s genome cannot be found" % genome)
        exit(1)
        
    aln_cmd = ["hisat2", "-x", genome, "--no-unal", "--no-hd", "-r", "-"]

    readsrecovered, totalreads = 0, len(read_chunk)  
    align_proc = subprocess.run(aln_cmd,
                                  input='\n'.join(read_chunk),
                                  capture_output = True,
                                  encoding="ascii")

    align_proc = align_proc.stdout.split("\n")

    for line in align_proc:
        if not line:
            continue
        readsrecovered += 1

        line = line.strip()
        cols = line.split()

        id_, flag, chrom, left, _, cigar_str = cols[:6]
        left  = int(left) - 1
        right = left

        mptr   = 0
        cigars = cigar_re.findall(cigar_str)
        cigars = [[cigar[-1], int(cigar[:-1])] for cigar in cigars]
        for cigar in cigars:
            cigar_op, leng = cigar
            if cigar_op in "M":
                right += leng
                for i in range(leng):
                    mptr += 1

            elif cigar_op in "ND":
                right += leng
                extend_meta(mptr, cigar_op, 1)

            elif cigar_op in "IS":
                right += leng
                for i in range(leng):
                    extend_meta(mptr, cigar_op, 1)
                    mptr += 1

            
            if cigar_op not in cigarcounts:
                cigarcounts[cigar_op] = 0

            cigarcounts[cigar_op] += 1

    print(metacigar) 
    print(cigarcounts)
    print(readsrecovered, totalreads)


def run_call(fnindex, files, out):
    sqidx = sqwho.SeqWho_Index(fnindex)
    sqidx.init_vectors(files)
    success = sqidx.run_vector_population()

    if not success:
        print("Error loading and running vectors. Exiting!")
        exit(1)

    sqidx.run_vector_test(out)

def main():
    parser = argparse.ArgumentParser(
        description='Seq-Who: Full Version Beta-1.0.0')
    parser.add_argument('-x', '--index',
                        dest     = 'idx',
                        type     = str,
                        required = True,
                        help     = 'Index name')
    parser.add_argument('-f', '--files',
                        dest     = 'files',
                        nargs    = '+',
                        required = True,
                        help     = 'Files to type')

    parser.add_argument('-o', '--out',
                        dest     = 'out',
                        type     = str,
                        default  = 'SeqWho_call',
                        help     = 'Name of output files (Default: "SeqWho_call"')

    args = parser.parse_args()

    ## Aligner options are not implemented yet
    # parser.add_argument('-a', '--aligner',
    #                     dest     = 'aligner',
    #                     default  = 'hisat2',
    #                     help     = 'Set the aligner to use \
    #                                 (Options: HISAT2, BWA, BOWTIE2)\
    #                                 (Default: HISAT2)')

    # args = parser.parse_args()
    # aligner = args.aligner.lower()
    # if not aligner in ['hisat2', 'bwa', 'bowtie2']:
    #     print('Error: selected aligner %s not supported' % aligner, 
    #           file=sys.stderr)
    #     exit(1)

    for f in args.files:
        if not os.path.exists(f):
            print("%s not found" % f)
            exit(0)

    run_call(args.idx, 
             args.files,
             args.out)

    return 0

if __name__ == '__main__':
    sys.exit(main())
