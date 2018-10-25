# nf-core/atacseq: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/). The initial QC and alignments are performed at the `library-level` e.g. if the same library has been sequenced more than once to increase sequencing depth. This has the advantage of being able to assess each library individually, and the ability to process multiple libraries from the same sample in parallel. The alignments are subsequently merged at the `replicate-level` and the `sample-level`. The latter involves merging the alignments across all replicates from the same experimental condition. This can be useful to increase the coverage for peak-calling and for other analyses that require high sequencing depth such as [motif footprinting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959825/).

* [Library-level analysis](#library-level-analysis)
    1. Raw read QC - [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    2. Adapter trimming - [`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
    3. Alignment - [`BWA`](https://sourceforge.net/projects/bio-bwa/files/)
    4. Mark duplicate reads - [`picard`](https://broadinstitute.github.io/picard/)
    5. Alignment filtering - [SAMTools](https://sourceforge.net/projects/samtools/files/samtools/), [BEDTools](https://github.com/arq5x/bedtools2/) & [Pysam](http://pysam.readthedocs.io/en/latest/installation.html

* [Replicate-level analysis](#replicate-level-analysis)
    1. Alignment merging - [`picard`](https://broadinstitute.github.io/picard/)
    2. Re-mark and remove duplicate reads - [`picard`](https://broadinstitute.github.io/picard/)
    3. Normalised bigWig files - [`BEDTools`](https://github.com/arq5x/bedtools2/) & [`wigToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/)
    4. TSS meta-profiles - [`deepTools`](https://deeptools.readthedocs.io/en/develop/)
    5. Call peaks - [`MACS2`](https://github.com/taoliu/MACS)
    6. Annotate peaks - [`HOMER`](http://homer.ucsd.edu/homer/download.html)
    7. Create consensus set of peaks - [`BEDTools`](https://github.com/arq5x/bedtools2/)
    8. Read counting relative to consensus set of peaks - [`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/)
    9. Differential binding analysis, PCA and clustering - [`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

> The analysis steps for the `sample-level` analysis are pretty much the same as for the `replicate-level` analysis. The only difference is that multiple libraries sequenced from the same sample will be merged at the `replicate-level` whereas all the replicates associated with an experimental condition will be merged at the `sample-level`.

* [Aggregate analysis](#aggregate-analysis)
    1. Collect and present QC at the raw read, alignment and peak-level - [`MultiQC`](http://multiqc.info/) & [`R`](https://www.r-project.org/)
    2. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation - [`IGV`](https://software.broadinstitute.org/software/igv/)
