#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on July 4th 2018 to create IGV session file from file list
#######################################################################
#######################################################################

import os
import argparse

import funcs

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Create a tab-delimited file with "file_path\tcolour" for IGV. This will be specific to directory structure of BABS-ATACSeqPE nextflow pipeline.'
Epilog = """Example usage: python igv_get_files.py <RESULTS_DIR> <OUT_FILE>"""
argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('RESULTS_DIR', help="Results directory. The directory structure used to find files will be specific to BABS-ATACSeqPE nextflow pipeline.")
argParser.add_argument('OUT_FILE', help="Path to output file.")
args = argParser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def igv_get_files(ResultsDir,OutFile):

    funcs.makedir(os.path.dirname(OutFile))

    ## GET SAMPLE-LEVEL FILES
    fileList = []
    for gid,odir,colList in [('mSm','sample',['255,0,0','0,0,178']),
                             ('mRp','replicate',['0,102,102','0,102,0'])]:

        fileList += [(x,'0,0,0') for x in funcs.recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), 'merged_peaks.%s.bed' % (gid))]
        fileList += [(x,colList[0]) for x in funcs.recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.FDR0.01.results.bed')]

        sampleFileDict = {}
        for ifile in funcs.recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.broadPeak') \
                   + funcs.recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.narrowPeak') \
                   + funcs.recursive_glob(os.path.join(ResultsDir,'bwa/','%s/' % (odir)), '*.bigWig'):
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
