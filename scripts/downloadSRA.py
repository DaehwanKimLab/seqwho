import pandas as pd
import os
from ftplib import FTP
import xml.etree.ElementTree as ET
import sys
from multiprocessing import Pool
import numpy as np
import string
from argparse import ArgumentParser
from itertools import repeat
import shutil
import subprocess
import time

########### Helper functions code from Seqmap/SRA_META/SRAParser.py ###########
class SampleParser:
    singleAttrib=['TITLE','DESCRIPTION','SCIENTIFIC_NAME','SUBMITTER_ID']
    attribName='SAMPLE_ATTRIBUTE'
    recordName='SAMPLE'
    f_extd='.sample.xml'
    def parseSample(self,sample):
        keys=[]
        vals=[]
        for myAttrib in self.singleAttrib:
            Iter=sample.iter(myAttrib)
            try:
                val=next(Iter).text
            except StopIteration:
                continue
            if val !=None and myAttrib !=None:
                keys.append(myAttrib),vals.append(val.replace('\t',''))
        for sample_attrib in sample.iter(self.attribName):
            val=sample_attrib.findtext('VALUE')
            if val!=None:
                keys.append(sample_attrib.findtext('TAG'))
                vals.append(val.replace('\t',''))
        return keys,vals
    def parseFile(self,fdir):
        cumAccessions,cumKeys,cumVals=[],[],[]
        root=ET.parse(fdir).getroot()
        for s in root.findall(self.recordName):
            accession=s.attrib['accession']
            keys,vals=self.parseSample(s)
            accessions=[accession]*len(keys)
            cumAccessions+=accessions
            cumKeys+=keys
            cumVals+=vals
        return cumAccessions,cumKeys,cumVals
    def parseFiles(self,sras,inDir):
        MultCumAccessions,MultCumKeys,MultCumVals=[],[],[]
        for sra in sras:
            fullDir=inDir+'/'+sra+'/'+sra+self.f_extd
            if os.path.isfile(fullDir):
                cumAccessions,cumKeys,cumVals=self.parseFile(fullDir)
                MultCumAccessions+=cumAccessions
                MultCumKeys+=cumKeys
                MultCumVals+=cumVals
        srsToKeyI=pd.MultiIndex.from_tuples(list(zip(MultCumAccessions,MultCumKeys)))
        multIS=pd.Series(data=MultCumVals,index=srsToKeyI)
        return multIS

class ExperimentParser:
    singleAttrib=['TITLE','DESIGN_DESCRIPTION','LIBRARY_STRATEGY','LIBRARY_SOURCE','LIBRARY_SELECTION']
    attribName='EXPERIMENT_ATTRIBUTE'
    recordName='EXPERIMENT'
    f_extd='.experiment.xml'
    def parseSample(self,sample):
        keys=[]
        vals=[]
        for myAttrib in self.singleAttrib:
            Iter=sample.iter(myAttrib)
            try:
                val=next(Iter).text
            except StopIteration:
                continue
            if val !=None and myAttrib !=None:
                keys.append(myAttrib),vals.append(val.replace('\t',''))
                #print myAttrib,val.replace('\t','')
        #extract the layout
        Iter=sample.iter('LIBRARY_LAYOUT')
        tmp=next(Iter)
        children=next(iter(tmp))
        keys.append('LIBRARY_LAYOUT'),vals.append(children.tag)
        #print children
        #print '\t'.join(map(str,[,children.text,children.items(),children.attrib,children.keys()]))
        #print
        for sample_attrib in sample.iter(self.attribName):
            val=sample_attrib.findtext('VALUE')
            if val!=None:
                keys.append(sample_attrib.findtext('TAG'))
                vals.append(val.replace('\t',''))
        return keys,vals
    def parseFile(self,fdir):
        cumAccessions,cumKeys,cumVals=[],[],[]
        root=ET.parse(fdir).getroot()
        for s in root.findall(self.recordName):
            accession=s.attrib['accession']
            keys,vals=self.parseSample(s)
            accessions=[accession]*len(keys)
            cumAccessions+=accessions
            cumKeys+=keys
            cumVals+=vals
        return cumAccessions,cumKeys,cumVals
    def parseFiles(self,sras,inDir):
        MultCumAccessions,MultCumKeys,MultCumVals=[],[],[]
        for sra in sras:
            fullDir=inDir+'/'+sra+'/'+sra+self.f_extd
            if os.path.isfile(fullDir):
                cumAccessions,cumKeys,cumVals=self.parseFile(fullDir)
                MultCumAccessions+=cumAccessions
                MultCumKeys+=cumKeys
                MultCumVals+=cumVals
        srsToKeyI=pd.MultiIndex.from_tuples(list(zip(MultCumAccessions,MultCumKeys)))
        multIS=pd.Series(data=MultCumVals,index=srsToKeyI)
        return multIS

class SRAParser:
    outputPostfix = '.srx.pickle' #just put one of them out there
    outputPostfix2='.srs.pickle'
    inputPostfix = ''
    untaredDir = ''

    def __init__(self):
        return


    def process(self,inDir):

        #global untaredDir
        untaredDir = '/'.join(inDir.split('/')[:-2]) + '/utaredDir/'
        nConcurJob = 2000

        chunkNum=int( inDir.split('/')[-1])#the number indicate the chunk of files to be processed

        allSra=os.listdir(untaredDir)
        chunkSize=int(len(allSra)/nConcurJob)
        #print (chunkSize)
        sras=allSra[(chunkNum*chunkSize):((chunkNum+1)*chunkSize)]
        #with open(inDir+self.outputPostfix,'w') as f:
        #print (sras)
        exprS = ExperimentParser().parseFiles(sras, untaredDir)
        sampleS=SampleParser().parseFiles(sras,untaredDir)

        exprS.to_pickle(inDir+self.outputPostfix)
        sampleS.to_pickle(inDir+self.outputPostfix2)


########## download SRA data ###########
# code from Seqmap/Pipelines/Update_SRA_meta_data/pull_SRA_meta.ipynb

def downloadSRAMeta(tmpDir):

    bdir = tmpDir
    #os.chdir(bdir)

    untaredDir = tmpDir + '/utaredDir/'
    if not os.path.exists(untaredDir):
        os.makedirs(untaredDir)

    ftpLink = 'ftp.ncbi.nlm.nih.gov'
    myRemoteDir = 'sra/reports/Metadata/'
    ftp = FTP(ftpLink)
    ftp.login()
    ftp.cwd(myRemoteDir)
    fnames = pd.Series(ftp.nlst())

    # sort the data numerically
    myFullSraMeta = fnames[fnames.str.contains('NCBI_SRA_Metadata_Full')].sort_values().iloc[-1]
    myDownloadFnames = ['SRA_Accessions.tab', myFullSraMeta,
                        'SRA_Run_Members.tab']

    print('Files to be downloaded from NCBI:', myDownloadFnames)

    for f in myDownloadFnames:
        fileDir = bdir + f
        File = open(fileDir, 'wb')
        ###reopen ftp everytime to avoid idling
        ftp = FTP(ftpLink)
        ftp.login()
        ftp.cwd(myRemoteDir)
        ftp.retrbinary('RETR %s' % f, File.write)
        File.close()

    tarCmd = 'tar --skip-old-files -xf {inDir} -C {out_dir}'.format(inDir=bdir + myFullSraMeta,
                                                                          out_dir=untaredDir)
    os.system(tarCmd)


########### Parse SRA data ####################
# code changed from Seqmap/SRA_META/SRAManager.ipynb

def prase(fname):
    p = SRAParser()
    p.process(fname)


def parseSRAData(parseDir, nConcurJob, nthreads):

    pool = Pool(nthreads)
    pool.map(prase, [parseDir + s for s in list(map(str, list(range(nConcurJob+1))))])


########## Merge SRA data ##########
# code from Seqmap/SRA_META/SRAmerge.ipynb

def mergeSRA(tmpDir, parseDir):

    global srsMergedS
    global srxMergedS

    #baseOutDir = tmpDir
    inDir = parseDir

    subsetIds = set(map(lambda s: s.split('.')[0], os.listdir(inDir)))

    srsS = []
    srxS = []
    for Id in subsetIds:
        srsS.append(pd.read_pickle(inDir + Id + '.srs.pickle'))
        srxS.append(pd.read_pickle(inDir + Id + '.srx.pickle'))
    srsMergedS = pd.concat(srsS)
    srxMergedS = pd.concat(srxS)
    #srsMergedS.to_pickle(baseOutDir + 'allSRS.pickle.gz')
    #srxMergedS.to_pickle(baseOutDir + 'allSRX.pickle.gz')


########## annotate SRA meta data based on SRX parsed ##########
# code from Seqmap/Pipelines/Update_SRA_meta_data/annotate_SRA_meta_data.ipynb

def annotateSRA(tmpDir, outputDir):

    global srsMergedS
    global srxMergedS

    basePickleDir = tmpDir
    accessionProjDir = tmpDir + '/SRA_Accessions.tab'
    accessionDir = tmpDir + '/SRA_Run_Members.tab'

    #srsMergedS = pd.read_pickle(basePickleDir + 'allSRS.pickle.gz')
    #srxMergedS = pd.read_pickle(basePickleDir + 'allSRX.pickle.gz')

    sra_dump_csv_dir = tmpDir + '/sra_dump.csv.gz'
    #sra_dump_pickle_dir = tmpDir + '/sra_dump.pickle'

    projAccessDf = pd.read_csv(accessionProjDir, sep='\t', dtype='str', index_col=0)

    accessionDf = pd.read_csv(accessionDir, sep='\t', dtype='str', index_col=0)

    validScientificNameSrs = srsMergedS[srsMergedS.index.get_level_values(1) == 'SCIENTIFIC_NAME']

    srsToSpecies = pd.Series(index=validScientificNameSrs.index.get_level_values(0), data=validScientificNameSrs.values)

    del srsMergedS

    uniqueSrsToSpeciesS = srsToSpecies.groupby(srsToSpecies.index).first()

    layoutS = srxMergedS[srxMergedS.index.get_level_values(1) == 'LIBRARY_LAYOUT'].groupby(level=0).first()

    libraryS = srxMergedS[srxMergedS.index.get_level_values(1) == 'LIBRARY_STRATEGY'].groupby(level=0).first()

    srxToSeqS = libraryS.groupby(level=0).first()

    del srxMergedS


    #validSrx = srxToSeqS.index.get_level_values(0)

    subAccessionDf = accessionDf[
        accessionDf.Sample.isin(uniqueSrsToSpeciesS.index) & accessionDf.Experiment.isin(srxToSeqS.index)].copy()

    del accessionDf

    subAccessionDf.loc[:, 'ScientificName'] = uniqueSrsToSpeciesS.loc[subAccessionDf.Sample.values].values[:]

    subAccessionDf['Run'] = subAccessionDf.index

    publicRunDf = subAccessionDf

    algnedSrx = srxToSeqS.loc[publicRunDf.Experiment].values
    srrToSeqS = pd.Series(data=algnedSrx, index=publicRunDf.index)
    publicRunDf.loc[:, 'LibraryStrategy'] = srrToSeqS[:]



    publicRunDf['LibraryLayout'] = layoutS.loc[publicRunDf['Experiment'].values].values

    projectSubDf = projAccessDf.loc[publicRunDf.index]
    projectSubDf.columns = 'proj_accession_' + projectSubDf.columns

    del projAccessDf

    yearAnnotatedPublicDf = pd.concat([publicRunDf, projectSubDf], axis=1)

    del subAccessionDf
    del publicRunDf
    del projectSubDf

    yearAnnotatedPublicDf['ScientificName'] = yearAnnotatedPublicDf['ScientificName'].str.replace(' ', '_')

    printable = set(string.printable)

    removeNonAscii = lambda s: "".join(filter(lambda x: x in printable, s))

    subExportDf = yearAnnotatedPublicDf.loc[:, yearAnnotatedPublicDf.columns != 'Run']

    subExportDf['new_ScientificName'] = subExportDf['ScientificName']

    del yearAnnotatedPublicDf

    m_4 = subExportDf['new_ScientificName'].str.contains('human')
    # Homo_sapiens
    subExportDf.loc[m_4, 'new_ScientificName'] = 'Homo_sapiens'
    # m_4=subExportDf['ScientificName'].str.contains('metagenome')
    m_5 = subExportDf['new_ScientificName'].str.contains('mouse')
    # Homo_sapiens
    subExportDf.loc[m_5, 'new_ScientificName'] = 'Mus_musculus'
    ###for soil, to eliminate human contamination, align to human reference
    m_6 = subExportDf['new_ScientificName'].str.contains('^soil_metagenome#')
    subExportDf.loc[m_6, 'new_ScientificName'] = 'Homo_sapiens'

    m_7 = subExportDf['new_ScientificName'].str.contains('^Canis')
    subExportDf.loc[m_7, 'new_ScientificName'] = 'Canis_familiaris'

    for myCol in [u'ScientificName', u'LibraryStrategy',
                  u'LibraryLayout', 'new_ScientificName']:
        subExportDf[myCol] = subExportDf[myCol].apply(removeNonAscii)

    #subExportDf.to_csv(sra_dump_csv_dir, compression='gzip')

    pickleExportInDf = subExportDf.copy()

    del subExportDf

    timeCols = ['proj_accession_Updated'
        , 'proj_accession_Published'
        , 'proj_accession_Received']
    for timeCol in timeCols:
        tmp_dateS = pickleExportInDf[timeCol].str.split('T').str[0]
        pickleExportInDf[timeCol] = pd.to_datetime(tmp_dateS, errors='ignore')

    ignoreCols = ['proj_accession_BioProject', 'proj_accession_Study', 'proj_accession_Spots',
                  'proj_accession_BioSample', 'proj_accession_Sample', 'proj_accession_Experiment',
                  'proj_accession_Alias', 'proj_accession_Md5sum',
                  'BioSample', 'proj_accession_Bases']

    subsetPickleDf = pickleExportInDf.loc[:, ~pickleExportInDf.columns.isin(ignoreCols)].copy()

    intTypes = ['Spots', 'Bases']
    for intType in intTypes:
        subsetPickleDf.loc[:, intType] = pd.to_numeric(subsetPickleDf.loc[:, intType], errors='coerce')

    factorTypes = [u'proj_accession_ReplacedBy', u'proj_accession_Type',
                   u'proj_accession_Status', u'Status', u'proj_accession_Visibility',
                   u'LibraryLayout', u'proj_accession_Loaded', u'LibraryStrategy',
                   u'proj_accession_Updated', u'proj_accession_Received',
                   u'proj_accession_Published', u'Member_Name', u'proj_accession_Center',
                   u'ScientificName', 'new_ScientificName']  # u'Study', u'proj_accession_Submission'

    for factorType in factorTypes:
        subsetPickleDf.loc[:, factorType] = subsetPickleDf.loc[:, factorType].astype('category')

    subsetPickleDf['Run_db'] = subsetPickleDf.index.str.extract('([EDS]RR)', expand=False).values
    subsetPickleDf['Run_digits'] = subsetPickleDf.index.str.extract('[EDS]RR(\d+)', expand=False).astype(np.int)

    subsetPickleDf = subsetPickleDf.sort_values(by=['ScientificName'])
    col = subsetPickleDf['ScientificName']
    ScientificNmaeDic = {}
    for count, name in enumerate(col):
        line = count + 1
        if name in ScientificNmaeDic:
            ScientificNmaeDic[name][1] = line
        else:
            ScientificNmaeDic[name] = [line, line]

    output_sra_dump_dir_index = outputDir + '/sra_dump_index.tsv'
    f_index = open(output_sra_dump_dir_index, 'w')
    for name, line in ScientificNmaeDic.items():
        f_index.write('%s\t%s\t%s\n' % (name, line[0], line[1]))
    f_index.close()

    output_sra_dump_dir = outputDir + '/sra_dump.tsv'
    subsetPickleDf.to_csv(output_sra_dump_dir, sep = '\t')

def removeDirIn(path):
    path += '/'
    for x in os.listdir(path):
        if os.path.isfile(path + x):
            continue
        cmd = "rm -rf %s" % (path + x)
        subprocess.Popen(cmd, shell=True)

def removeFileIn(path):
    path += '/'
    for x in os.listdir(path):
        if os.path.isfile(path + x):
            cmd = "rm %s" % (path + x)
            subprocess.Popen(cmd, shell=True)

def checkDirExistIn(path):
    # return True is there is dir in path.
    # return False is there is no dir in path.
    path += '/'
    for x in os.listdir(path):
        if not os.path.isfile(path + x):
            return True
    return False


def runAll(outputDir, nthreads):

    # untaredDir = tmpDir + '/utaredDir/'
    tmpDir = outputDir + 'tmp_SRA/'
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)

    parseDir = tmpDir + '/SRA_parse/'
    if not os.path.exists(parseDir):
        os.makedirs(parseDir)

    nConcurJob = 2000

    print("downloading SRA data")
    downloadSRAMeta(tmpDir)

    print("parsing SRA data")
    parseSRAData(parseDir, nConcurJob, nthreads)

    print("merging SRA data")
    mergeSRA(tmpDir, parseDir)

    removeDirIn(tmpDir)

    print("annotating SRA data")
    annotateSRA(tmpDir, outputDir)

    removeFileIn(tmpDir)

    while (len(os.listdir(tmpDir)) != 0):
        time.sleep(5)

    shutil.rmtree(tmpDir)


if __name__ == '__main__':
    parser = ArgumentParser(description='Download SRA database')

    parser.add_argument("--out-dir",
                        "-o",
                        dest="out_dir",
                        type=str,
                        default=".")
    parser.add_argument("--threads",
                        "-p",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    args = parser.parse_args()

    outputDir = args.out_dir + '/'
    nthreads = args.threads

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    srsMergedS = None
    srxMergedS = None

    runAll(outputDir, nthreads)



