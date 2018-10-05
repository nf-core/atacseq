#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on July 4th 2018 to create IGV session file from file list
#######################################################################
#######################################################################

import os
import argparse
import subprocess
import fnmatch

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Create a tab-delimited file with "file_path\tcolour" for IGV. This will be specific to results directory structure of nf-core/atacseq nextflow pipeline.'
Epilog = """Example usage: python igv_get_files.py <RESULTS_DIR> <OUT_FILE>"""
argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('RESULTS_DIR', help="Results directory. The results directory structure used to find files will be specific to nf-core/atacseq nextflow pipeline.")
argParser.add_argument('OUT_FILE', help="Path to output file.")
args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
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
############################################
## MAIN FUNCTION
############################################
############################################

def igv_get_files(ResultsDir,OutFile):

    makedir(os.path.dirname(OutFile))

    ## GET SAMPLE-LEVEL FILES
    fileList = []
    for gid,odir,colList in [('mSm','sample',['255,0,0','0,0,178']),
                             ('mRp','replicate',['0,102,102','0,102,0'])]:

        fileList += [(x,'0,0,0') for x in recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), 'merged_peaks.%s.bed' % (gid))]
        fileList += [(x,colList[0]) for x in recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.FDR0.01.results.bed')]

        sampleFileDict = {}
        for ifile in recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.broadPeak') \
                   + recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.narrowPeak') \
                   + recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.bigWig'):
            extension = os.path.splitext(ifile)[1]
            sampleid = ''
            if extension == '.broadPeak':
                sampleid = os.path.basename(ifile).replace('_peaks.broadPeak','')
            elif extension == '.narrowPeak':
                sampleid = os.path.basename(ifile).replace('_peaks.narrowPeak','')
            elif extension == '.bigWig':
                sampleid = os.path.basename(ifile).replace('.bigWig','')
            if not sampleFileDict.has_key(sampleid):
                sampleFileDict[sampleid] = []
            sampleFileDict[sampleid].append((ifile,colList[1]))
        for sampleid in sorted(sampleFileDict.keys()):
            fileList += sampleFileDict[sampleid]

    fout = open(OutFile,'w')
    for ifile,colour in fileList:
        fout.write('%s\t%s\n' % (ifile,colour))
    fout.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

igv_get_files(ResultsDir=args.RESULTS_DIR,OutFile=args.OUT_FILE)

############################################
############################################
############################################
############################################
