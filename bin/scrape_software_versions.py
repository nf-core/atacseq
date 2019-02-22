#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/atacseq': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'BamTools': ['v_bamtools.txt', r"bamtools (\S+)"],
    'Picard': ['v_picard.txt', r"([\d\.]+)-SNAPSHOT"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'Pysam': ['v_pysam.txt', r"(\S+)"],
    'MACS2': ['v_macs2.txt', r"macs2 (\S+)"],
    'HOMER': ['v_homer.txt', r"(\S+)"],
    'ataqv': ['v_ataqv.txt', r"(\S+)"],
    'featureCounts': ['v_featurecounts.txt', r"featureCounts v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}

results = OrderedDict()
results['nf-core/atacseq'] = 'NA'
results['Nextflow'] = 'NA'
results['FastQC'] = 'NA'
results['Trim Galore!'] = 'NA'
results['BWA'] = 'NA'
results['Samtools'] = 'NA'
results['BEDTools'] = 'NA'
results['BamTools'] = 'NA'
results['Picard'] = 'NA'
results['R'] = 'NA'
results['Pysam'] = 'NA'
results['MACS2'] = 'NA'
results['HOMER'] = 'NA'
results['ataqv'] = 'NA'
results['featureCounts'] = 'NA'
results['MultiQC'] = 'NA'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to TSV
for k,v in results.items():
    print("{}\t{}".format(k,v))
