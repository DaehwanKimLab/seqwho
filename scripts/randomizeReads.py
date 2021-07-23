import sys
import os
import random
import gzip
import io

def readFastq(fname):
    if not os.path.exists(fname):
        print("No file with name %s" % fname)
        exit(0)

    fullFastq = []
    ixFastq   = []
    fastqLine = []
    counter   = 0
    with gzip.open(fname, 'rb') as ifi:
        for line in ifi:
            if line.decode("utf-8").startswith("@"):
                if fastqLine:
                    ixFastq.append(counter)
                    fullFastq.append(fastqLine)
                    counter += 1

                fastqLine = []
            fastqLine.append(line.strip())
        ixFastq.append(counter)
        fullFastq.append(fastqLine)
    return ixFastq, fullFastq

def writeFastq(fname, ixFastq, fullFastq):
    os.remove(fname)
    with gzip.open(fname, 'wb') as ofo:
        for ix in ixFastq:
            for line in fullFastq[ix]:
                ofo.write(line)
                ofo.write('\n'.encode("utf-8"))

def randomizeFastq(seed, fname):
    try:
        ix, fastq = readFastq(fname)
        random.seed(seed)
        random.shuffle(ix)
        print(ix[0:10])
        writeFastq(fname, ix, fastq)
    except:
        print("Failed to randomize %s" % fname)
        exit(0)

if __name__ == '__main__':
    fname = sys.argv[1]
    seed  = sys.argv[2]
    randomizeFastq(seed, fname)