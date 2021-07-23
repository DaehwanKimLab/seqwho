import pandas as pd
from argparse import ArgumentParser


def filterByDate(data, argument):
    data = data[data['proj_accession_Published'] != '-']

    if argument[0] == '=':
        date = argument[1:]
        data = data[data['proj_accession_Published'] == date]
        return data
    if argument[:2] == '>=':
        date = argument[2:]
        data = data[data['proj_accession_Published'] >= date]
        return data
    if argument[:2] == '<=':
        date = argument[2:]
        data = data[data['proj_accession_Published'] <= date]
        return data
    if argument[0] == '>':
        date = argument[1:]
        data = data[data['proj_accession_Published'] > date]
        return data
    if argument[0] == '<':
        date = argument[1:]
        data = data[data['proj_accession_Published'] < date]
        return data


def filterByColumn(data, argument, columnName):
    argumentList = argument.split(',')
    filterList = False * len(data)

    while (len(argumentList) != 0):
        argumentValue = argumentList[0]
        argumentList = argumentList[1:]
        filterList = (data[columnName] == argumentValue) | filterList

    data = data[filterList]

    return data


def filterByIndex(data, argument):
    argumentList = argument.split(',')

    data = data.loc[argumentList]

    return data


def filterByValue(data, argument, columnName):
    #data = data[data[columnName] != '-']

    if 'AND' in argument:
        argumentList = argument.split('AND')
        filterList = True * len(data)

        while (len(argumentList) != 0):
            argumentValue = argumentList[0]
            argumentList = argumentList[1:]

            if argumentValue[0] == '=':
                value = int(argumentValue[1:])
                filterList = (data[columnName] == value) & filterList
                continue
            if argumentValue[:2] == '>=':
                value = int(argumentValue[2:])
                filterList = (data[columnName] >= value) & filterList
                continue
            if argumentValue[:2] == '<=':
                value = int(argumentValue[2:])
                filterList = (data[columnName] <= value) & filterList
                continue
            if argumentValue[0] == '>':
                value = int(argumentValue[1:])
                filterList = (data[columnName] > value) & filterList
                continue
            if argumentValue[0] == '<':
                value = int(argumentValue[1:])
                filterList = (data[columnName] > value) & filterList
                continue

    else:
        argumentList = argument.split('OR')
        filterList = False * len(data)

        while (len(argumentList) != 0):
            argumentValue = argumentList[0]
            argumentList = argumentList[1:]

            if argumentValue[0] == '=':
                value = int(argumentValue[1:])
                filterList = (data[columnName] == value) | filterList
                continue
            if argumentValue[:2] == '>=':
                value = int(argumentValue[2:])
                filterList = (data[columnName] >= value) | filterList
                continue
            if argumentValue[:2] == '<=':
                value = int(argumentValue[2:])
                filterList = (data[columnName] <= value) | filterList
                continue
            if argumentValue[0] == '>':
                value = int(argumentValue[1:])
                filterList = (data[columnName] > value) | filterList
                continue
            if argumentValue[0] == '<':
                value = int(argumentValue[1:])
                filterList = (data[columnName] > value) | filterList
                continue

    data = data[filterList]
    return data

def filterSRA(data, species, experiment, sample, run, study, spots, bases, date, all_, layout, seq):

    if not all_:
        data = data[data['proj_accession_Visibility'] == 'public']

    if species:
        data = filterByColumn(data, species, 'ScientificName')

    if experiment:
        data = filterByColumn(data, experiment, 'Experiment')

    if sample:
        data = filterByColumn(data, sample, 'Sample')

    if run:
        data = filterByIndex(data, run)

    if spots:
        data = filterByValue(data, spots, 'Spots')

    if bases:
        data = filterByValue(data, bases, 'Bases')

    if date:
        data = filterByDate(data, date)

    if study:
        data = filterByColumn(data, study, 'Study')

    if seq:
        data = filterByColumn(data, seq, 'LibraryStrategy')

    if layout:
        if layout == 'single':
            data = filterByColumn(data, 'SINGLE', 'LibraryLayout')
        elif layout == 'double':
            data = filterByColumn(data, 'PAIRED', 'LibraryLayout')
        else:
            print("Please enter 'single' or 'double' for --layout argument."
                  "The data is not filtered by layout argument for this time.")
    return data


if __name__ == '__main__':
    parser = ArgumentParser(description='SRA database Filter')

    parser.add_argument("--in-dir",
                        "-i",
                        dest="in_dir",
                        type=str,
                        help="Please enter the directory for the SRA database tsv file.")

    parser.add_argument("--out-dir",
                        "-o",
                        dest="out_dir",
                        type=str,
                        default=None,
                        help="Please enter the output directory include file name. If no argument entered."
                             "the output file name will be filtered_SRA.tsv under input directory.")

    parser.add_argument("--species",
                        "-sp",
                        dest="species",
                        type=str,
                        default=None,
                        help="If you enter more than one species, please seperate them by comma.\
                        please use _ to seperate the word for each species\
                        for example: -sp Homo_sapiens,Passerella_iliaca")

    parser.add_argument("--experiment",
                        "-e",
                        dest="experiment",
                        type=str,
                        default=None,
                        help="If you enter more than one experiment accession numbers,\
                        please seperate them by comma.\
                        for example: -e SRX3204478,SRX3204539")

    parser.add_argument("--sample",
                        "-sa",
                        dest="sample",
                        type=str,
                        default=None,
                        help="If you enter more than one sample accession numbers,\
                        please seperate them by comma.\
                        for example: -sa SRS2530621,SRS3932940")

    parser.add_argument("--run",
                        "-r",
                        dest="run",
                        type=str,
                        default=None,
                        help="If you enter more than one run accession numbers,\
                        please seperate them by comma.\
                        for example: -r ERR925021,SRR9611265")

    parser.add_argument("--spots",
                        "-spots",
                        dest="spots",
                        type=str,
                        default=None,
                        help="You can choose to filter the number of spots by =, >, >=, <, <=. \
                        if you have more than one condition, please join them by AND or OR.\
                        Please use single quotation for your argument, for example: -spot '>100000OR<=800000' \
                        Please do NOT put AND and OR in same argument")

    parser.add_argument("--bases",
                        "-b",
                        dest="bases",
                        type=str,
                        default=None,
                        help="You can choose to filter the number of bases by =, >, >=, <, <=. \
                        if you have more than one condition, please join them by AND or OR.\
                        Please use single quotation for your argument, for example: -b '>100000OR<=800000' \
                        Please do NOT put AND and OR in same argument")

    parser.add_argument("--date",
                        "-d",
                        dest="date",
                        type=str,
                        default=None,
                        help="You can choose to filter the publication date by =, >, >=, <, <=. \
                        Please use YYYY-MM-DD as date format. You may ignore DD or MM-DD. \
                        Please use single quotation for your argument, for example: \
                        enter -d '>2018-09' to get the SRA published after Sep/2018. \
                        enter -d '<= 2017-07-13' to get the SRA published before or on 7/13/2017.\
                        enter -d '=2016-07-03' to get the SRA published on 7/3/2016")

    parser.add_argument("--study",
                        "-st",
                        dest="study",
                        type=str,
                        default=None,
                        help="If you enter more than one study accession numbers,\
                        please seperate them by comma.\
                        for example: -st SRP212249,SRP212314")

    parser.add_argument("--all",
                        "-all",
                        dest="all_",
                        action="store_true",
                        help="Output both controlled access and public accession. Output will only\
                        contain public accession by defualt")

    parser.add_argument("--layout",
                        "-la",
                        dest="layout",
                        type=str,
                        default=None,
                        help="You can choose whether you want single end or paired end data."
                             "Please enter -l single for single end data only,"
                             "-l double for paired end data only. If there is no input for this argument,"
                             "output data will include both single or paired end data")

    parser.add_argument("--seq",
                        "-seq",
                        dest="seq",
                        type=str,
                        default=None,
                        help="Please choose Sequencing stategy from: \
                        AMPLICON,ATAC-seq,Bisulfite-Seq,CLONE,CLONEEND,CTS,ChIA-PET,ChIP,\
                        ChIP-Seq,DNase-Hypersensitivity,EST,FAIRE-seq,FINISHING,FL-cDNA,Hi-C,\
                        MBD-Seq,MNase-Seq,MRE-Seq,MeDIP-Seq,OTHER,POOLCLONE,RAD-Seq,RIP-Seq,\
                        RNA-Seq,SELEX,Synthetic-Long-Read,Targeted-Capture,\
                        Tethered Chromatin Conformation Capture,Tn-Seq,VALIDATION,\
                        WCS,WGA,WGS,WXS,miRNA-Seq,ncRNA-Seq,other.  \
                        Please seperate them by comma and do NOT put space in between.\
                        For example: -seq ChIP-Seq,WGA,WGS")

    args = parser.parse_args()

    fileName = args.in_dir
    outputDir = args.out_dir
    species = args.species
    experiment = args.experiment
    sample = args.sample
    run = args.run
    spots = args.spots
    bases = args.bases
    date = args.date
    study = args.study
    all_ = args.all_
    seq = args.seq
    layout = args.layout

    data = pd.read_csv(fileName, sep='\t', index_col=0)

    filteredData = filterSRA(data, species, experiment, sample, run, study, spots, bases, date, all_, layout, seq)

    if not outputDir:
        outputDir = fileName.split('/')
        outputDir[-1] = 'filtered_SRA.tsv'
        outputDir = '/'.join(outputDir)

    filteredData.to_csv(outputDir, sep='\t')