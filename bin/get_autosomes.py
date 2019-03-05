#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on February 26th 2019 to render methods
#######################################################################
#######################################################################

import os
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Get a list of autosomes from FAI file for assembly and write to file.'
Epilog = """Example usage: python get_autosomes.py <FAI_FILE> <OUT_FILE>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('METHODS', help="FAI input file.")
argParser.add_argument('OUT_FILE', help="Output file containing one chromosome per line.")
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
############################################
## MAIN FUNCTION
############################################
############################################

## TESTED WITH IGENOMES ASSEMBLIES FOR:
    # Arabidopsis_thaliana
    # Caenorhabditis_elegans
    # Danio_rerio
    # Drosophila_melanogaster
    # Enterobacteriophage_lambda
    # Escherichia_coli_K_12_DH10B
    # Escherichia_coli_K_12_MG1655
    # Gallus_gallus
    # Homo_sapiens
    # Mus_musculus
    # Mycobacterium_tuberculosis_H37RV
    # PhiX
    # Rattus_norvegicus
    # Saccharomyces_cerevisiae
    # Schizosaccharomyces_pombe
    # Staphylococcus_aureus_NCTC_8325

def get_autosomes(FAIFile,OutFile):

    makedir(os.path.dirname(OutFile))

    ## READ IN CHROMOSOME IDS
    chrList = []
    fin = open(FAIFile,'r')
    while True:
        line = fin.readline()
        if line:
            chrList.append(line.strip().split('\t')[0])
        else:
            fin.close()
            break

    ## REMOVE EXACT MATCHES TO MITOCHONDRIAL, SEX AND OTHER NON-AUTOSOMAL CHROMOSOMES
    exactMatch = ['y', 'w', 'z', 'm', 'mt', 'mtdna', 'pt', 'mtr', '2-micron', 'ebv']
    if 'VI' not in chrList and 'chrVI' not in chrList:                  ## Caenorhabditis elegans has roman numerals and chrX but only 5 chromosomes!
        exactMatch += ['x']
    filteredList = [x for x in chrList if x.lower() not in exactMatch and x.lower().lstrip('chr') not in exactMatch]

    ## REMOVE INEXACT MATCHS TO SMALLER/RANDOM HUMAN CONTIGS
    fuzzyMatch = ['chrUn', 'random', 'v1', 'v2', 'kn707', 'jtfh']
    filteredList = [x for x in filteredList if not any(y in x.lower() for y in fuzzyMatch)]

    ## WRITE TO FILE
    fout = open(OutFile,'w')
    for chrom in filteredList:
        fout.write('%s\n' % (chrom))
    fout.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

get_autosomes(FAIFile=args.FAI_FILE,OutFile=args.OUT_FILE)

############################################
############################################
############################################
############################################
