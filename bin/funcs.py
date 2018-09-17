#######################################################################
#######################################################################
## Created on July 6th 2018 to store utility functions for .py exe
#######################################################################
#######################################################################

import os
import fnmatch
import subprocess

############################################
############################################
## GENERAL FUNCTIONS
############################################
############################################

def makedir(path):

    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################

def recursive_glob(treeroot, pattern):

    results = []
    for base, dirs, files in os.walk(os.path.abspath(treeroot)):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)

    return results

############################################

def numLinesInFile(File):

    p = subprocess.Popen(['wc', '-l', File], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

############################################

def intersect(a,b):

    return filter(lambda x:x in a,b)

############################################

def percentToStr(numer,denom,sigFigs=2,parentheses=False):

    percStr = ''
    sigFigStr = '%.' + str(sigFigs) + 'f'
    if numer == 0 or denom == 0:
        percStr = sigFigStr % (0.0) + '%'
    elif numer == 0 and denom == 0:
        percStr = sigFigStr % (0.0) + '%'
    else:
        percStr = sigFigStr %((numer/float(denom))*100) + '%'
    if parentheses != False:
        return '(' + percStr + ')'
    else:
        return percStr

############################################
############################################
## SOFTWARE FUNCTIONS
############################################
############################################

def cutadaptPELogToDict(cutadaptLogFile):

    cutadaptDict = {}
    fin = open(cutadaptLogFile,'r')
    fields = ['Total read pairs processed:','','','']
    for line in fin.readlines():
        if line.find('Total read pairs processed:') != -1:
            cutadaptDict['totalPairs'] = int(line.split(':')[1].strip().replace(',',''))
        elif line.find('Pairs written (passing filters):') != -1:
            cutadaptDict['passTrimmedPairs'] = int(line.split(':')[1].strip().split()[0].replace(',',''))
        elif line.find('Total written (filtered):') != -1:
            cutadaptDict['passTrimmedBases'] = line.split(':')[1].split()[-1][1:-1]
    fin.close()

    return cutadaptDict

############################################

def flagstatToDict(flagStatFile):

    flagstatDict = {}
    fin = open(flagStatFile,'r')
    for line in fin.readlines():
        lspl = line.split('(')[0].split()
        flagstatDict[' '.join(lspl[3:])] = (int(lspl[0]),int(lspl[2]))
    fin.close()

    return flagstatDict

############################################

def idxstatsToDict(idxstatsFile):

    idxstatsDict = {}
    fin = open(idxstatsFile,'r')
    for line in fin.readlines():
        chrom,clen,mapped,unmapped = line.strip().split('\t')
        idxstatsDict[chrom] = (int(clen),int(mapped),int(unmapped))
    fin.close()

    return idxstatsDict

############################################

def picardInsertMetricsToDict(insertMetricsFile):

    metricsDict = {}
    fin = open(insertMetricsFile,'r')
    lines = fin.readlines()
    for idx in range(len(lines)):
        if lines[idx][:len('MEDIAN_INSERT_SIZE')] == 'MEDIAN_INSERT_SIZE':
            metricsDict = dict(zip(lines[idx].strip().split('\t'),lines[idx+1].strip().split('\t')))
            fin.close()
            break

    return metricsDict

############################################
############################################
############################################
############################################
