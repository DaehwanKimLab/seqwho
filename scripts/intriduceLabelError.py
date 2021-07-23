import sys
import os
import random

def load_labels(fname):
    if(not os.path.exists(fname)):
        print("No file found with %s" % fname)
        exit(1)

    fileLines = open(fname, "r").read()
    fileLines = fileLines.split("\n")

    processedFileLines = {"Lines"   : [],
                          "Species" : [],
                          "Library" : []}
    for i in range(len(fileLines)):
        fL = fileLines[i].split(",")
        if len(fL) != 3:
            continue
        if fL[1] not in processedFileLines["Species"]:
            processedFileLines["Species"].append(fL[1])
        if fL[2] not in processedFileLines["Library"]:
            processedFileLines["Library"].append(fL[2])

        processedFileLines["Lines"].append(fL)

    processedFileLines["Species"] = sorted(processedFileLines["Species"])
    processedFileLines["Library"] = sorted(processedFileLines["Library"])

    print("%d entried were loaded and processed" % len(processedFileLines["Lines"]))
    return processedFileLines

def printErrorFile(fname, lines):
    with open(fname, 'w') as ofo:
        for line in lines:
            ofo.write(",".join(line))
            ofo.write("\n")

def addErrorMixedSpecLib(fname, seedbase, mixmode):
    validmixmode = {"all", "species", "library"}
    if mixmode not in validmixmode:
        print("Invalid Mix mode. Use following options: %s" % str(validmixmode))
        exit(0)

    fileLines = load_labels(fname)
    maxVal    = len(fileLines["Lines"])
    totalDraw = int(maxVal*0.05) # 5% draw from fileLines lines

    random.seed(seedbase)
    fLindex = random.sample(range(0, maxVal), totalDraw)

    random.seed(seedbase)
    if mixmode == 'all':
        fLbool = [random.randint(0, 1) for i in range(totalDraw)]
    elif mixmode == 'species':
        fLbool = [1 for i in range(totalDraw)]
    else:
        fLbool = [0 for i in range(totalDraw)]

    setDraw = int(totalDraw/5)
    cnt    = 0
    cycles = 1
    for i in range(len(fLindex)):
        ix   = fLindex[i]
        spec = False
        if fLbool[i] == 1:
            spec = True

        print(fileLines["Lines"][ix])
        if spec:
            # Swap Species and assumes only two species
            species = fileLines["Lines"][ix][1]
            for j in range(len(fileLines["Species"])):
                if species != fileLines["Species"][j]:
                    fileLines["Lines"][ix][1] = fileLines["Species"][j]
        else:
            # Randomly select alternative library type
            library  = fileLines["Lines"][ix][2]
            selection = [x for x in fileLines["Library"] if library != x] # Make new selection list without the current library selection
            random.seed(seedbase)
            selIx = random.randint(0, len(selection) - 1)
            fileLines["Lines"][ix][2] = selection[selIx]
            seedbase += 1

        print(fileLines["Lines"][ix])
        print("<<<<<")
        cnt += 1
        if cnt == setDraw:
            printErrorFile(fname.replace(".txt", "_" + str(cycles) + "_" + mixmode + ".txt"), fileLines["Lines"])
            cnt    = 0
            cycles += 1
            


if __name__ == '__main__':
    fname    = "/Users/chrisbennett/OneDrive - University of Texas Southwestern/Publications/SeqWho_bioinformatics_202106_majorRevisions/extra_data/sraSeqWhoLabels16k.txt"
    seedbase = 100
    mode     = sys.argv[1]
    #fname    = sys.argv[2]
    #seedbase = sys.argv[3]

    
    print(fname, seedbase, mode)
    addErrorMixedSpecLib(fname, seedbase, mode)
