# ![nfcore/atacseq](docs/images/nf-core-atacseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/atacseq.svg?branch=master)](https://travis-ci.org/nf-core/atacseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/atacseq.svg)](https://hub.docker.com/r/nfcore/atacseq)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2634133.svg)](https://doi.org/10.5281/zenodo.2634133)

## Introduction

**nfcore/atacseq** is a bioinformatics analysis pipeline used for ATAC-seq data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Filtering to remove:
        * reads mapping to mitochondrial DNA ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); *paired-end only*)
        * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
    3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
    4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
    6. Calculate genome-wide enrichment ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
    7. Call broad/narrow peaks ([`MACS2`](https://github.com/taoliu/MACS))
    8. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    9. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    10. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    11. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    12. Generate ATAC-seq specific QC html report ([`ataqv`](https://github.com/ParkerLab/ataqv))
6. Merge filtered alignments across replicates ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Remove duplicate reads ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    4. Call broad/narrow peaks ([`MACS2`](https://github.com/taoliu/MACS))
    5. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    6. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    7. Count reads in consensus peaks relative to merged library-level alignments ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    8. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
8. Present QC for raw read, alignment, peak-calling and differential accessibility results ([`ataqv`](https://github.com/ParkerLab/ataqv), [`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/atacseq -profile test,<docker/singularity/conda>
```

iv. Start running your own analysis!

```bash
nextflow run nf-core/atacseq -profile <docker/singularity/conda> --input design.csv --genome GRCh37
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation
The nf-core/atacseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

The pipeline was originally written by [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) for use at [The Francis Crick Institute](https://www.crick.ac.uk/), London.

The pipeline was developed by [Harshil Patel](mailto:harshil.patel@crick.ac.uk).

The [nf-core/rnaseq](https://github.com/nf-core/rnaseq) and [nf-core/chipseq](https://github.com/nf-core/chipseq) pipelines developed by Phil Ewels were initially used as a template for this pipeline. Many thanks to Phil for all of his help and advice, and the team at SciLifeLab.

Many thanks to others who have helped out along the way too, including (but not limited to): [@apeltzer](https://github.com/apeltzer), [@sven1103](https://github.com/sven1103), [@MaxUlysse](https://github.com/MaxUlysse), [@micans](https://github.com/micans), [@jinmingda](https://github.com/jinmingda), [@ktrns](https://github.com/ktrns), [@crickbabs](https://github.com/crickbabs), [@pditommaso](https://github.com/pditommaso).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/atacseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citation

If you use nf-core/atacseq for your analysis, please cite it using the following doi: [10.5281/zenodo.2634132](https://doi.org/10.5281/zenodo.2634132)

You can cite the `nf-core` pre-print as follows:  
> Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v1).

## Software Citation

### Pipeline tools

[BWA](https://www.ncbi.nlm.nih.gov/pubmed/19451168/)
> Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009 Jul 15;25(14):1754-60. doi: 10.1093/bioinformatics/btp324. Epub 2009 May 18. PubMed PMID: 19451168; PubMed Central PMCID: PMC2705234.

[BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/20110278/)
> Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PubMed PMID: 20110278; PubMed Central PMCID: PMC2832824.

[SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943/)
> Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PubMed PMID: 19505943; PubMed Central PMCID: PMC2723002.

[BamTools](https://www.ncbi.nlm.nih.gov/pubmed/21493652/)
> Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011 Jun 15;27(12):1691-2. doi: 10.1093/bioinformatics/btr174. Epub 2011 Apr 14. PubMed PMID: 21493652; PubMed Central PMCID: PMC3106182.

[UCSC tools](https://www.ncbi.nlm.nih.gov/pubmed/20639541/)
> Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics. 2010 Sep 1;26(17):2204-7. doi: 10.1093/bioinformatics/btq351. Epub 2010 Jul 17. PubMed PMID: 20639541; PubMed Central PMCID: PMC2922891.

[preseq](https://www.ncbi.nlm.nih.gov/pubmed/23435259/)
> Daley T, Smith AD. Predicting the molecular complexity of sequencing libraries. Nat Methods. 2013 Apr;10(4):325-7. doi: 10.1038/nmeth.2375. Epub 2013 Feb 24. PubMed PMID: 23435259; PubMed Central PMCID: PMC3612374.

[deepTools](https://www.ncbi.nlm.nih.gov/pubmed/27079975/)
> Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257. Epub 2016 Apr 13. PubMed PMID: 27079975; PubMed Central PMCID: PMC4987876.

[MACS2](https://www.ncbi.nlm.nih.gov/pubmed/18798982/)
> Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PubMed PMID: 18798982; PubMed Central PMCID: PMC2592715.

[HOMER](https://www.ncbi.nlm.nih.gov/pubmed/20513432/)
> Heinz S, Benner C, Spann N, Bertolino E, Lin YC, Laslo P, Cheng JX, Murre C, Singh H, Glass CK. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell. 2010 May 28;38(4):576-89. doi: 10.1016/j.molcel.2010.05.004. PubMed PMID: 20513432; PubMed Central PMCID: PMC2898526.

[featureCounts](https://www.ncbi.nlm.nih.gov/pubmed/24227677/)
> Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PubMed PMID: 24227677.

[DESeq2](https://www.ncbi.nlm.nih.gov/pubmed/25516281/)
> Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. PubMed PMID: 25516281; PubMed Central PMCID: PMC4302049.

[MultiQC](https://www.ncbi.nlm.nih.gov/pubmed/27312411/)
> Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[picard-tools](http://broadinstitute.github.io/picard)

[pysam](https://github.com/pysam-developers/pysam)

[ataqv](https://github.com/ParkerLab/ataqv)

### R packages

[R](https://www.R-project.org/)
> R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

[optparse](https://CRAN.R-project.org/package=optparse)
> Trevor L Davis (2018). optparse: Command Line Option Parser.

[RColorBrewer]: https://CRAN.R-project.org/package=RColorBrewer
> Erich Neuwirth (2014). RColorBrewer: ColorBrewer Palettes.

[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
> H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

[reshape2](http://www.jstatsoft.org/v21/i12/)
> Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20.

[scales](https://CRAN.R-project.org/package=scales)
> Hadley Wickham (2018). scales: Scale Functions for Visualization.

[pheatmap](https://CRAN.R-project.org/package=pheatmap)
> Raivo Kolde (2018). pheatmap: Pretty Heatmaps.

[lattice](https://cran.r-project.org/web/packages/lattice/index.html)
> Sarkar, Deepayan (2008) Lattice: Multivariate Data Visualization with R. Springer, New York. ISBN 978-0-387-75968-5

[vsn](https://bioconductor.org/packages/release/bioc/html/vsn.html)
> Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie Poustka and Martin Vingron. Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression. Bioinformatics 18, S96-S104 (2002).

[UpSetR](https://CRAN.R-project.org/package=UpSetR)
> Nils Gehlenborg (2017). UpSetR: A More Scalable Alternative to Venn and Euler Diagrams for Visualizing Intersecting Sets.

[xfun](https://CRAN.R-project.org/package=xfun)
> Yihui Xie (2018). xfun: Miscellaneous Functions by 'Yihui Xie'.

### Infrastructure tools

[Bioconda](https://www.ncbi.nlm.nih.gov/pubmed/29967506/)
> Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

[Anaconda](https://anaconda.com)
> Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

[Singularity](https://www.ncbi.nlm.nih.gov/pubmed/28494014/)
> Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
