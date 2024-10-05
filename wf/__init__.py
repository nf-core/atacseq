from typing import List, Optional

from latch import map_task
from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

import wf.Post_Process_Tasks as pp
from wf.entrypoint import Aligner, Reference, SampleSheet, initialize, nextflow_runtime


@workflow(metadata._nextflow_metadata)
def nf_nf_core_atacseq(
    input: List[SampleSheet],
    run_name: str,
    seq_center: Optional[str],
    read_length: Optional[int],
    with_control: bool,
    email: Optional[str],
    multiqc_title: Optional[str],
    genome_source: str,
    genome: Optional[str],
    latch_genome: Optional[Reference],
    fasta: Optional[LatchFile],
    gtf: Optional[LatchFile],
    gff: Optional[LatchFile],
    bwa_index: Optional[str],
    bowtie2_index: Optional[LatchDir],
    chromap_index: Optional[str],
    star_index: Optional[str],
    gene_bed: Optional[LatchFile],
    tss_bed: Optional[LatchFile],
    macs_gsize: Optional[float],
    blacklist: Optional[str],
    mito_name: Optional[str],
    save_reference: bool,
    ataqv_mito_reference: Optional[str],
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    skip_trimming: bool,
    save_trimmed: bool,
    keep_dups: bool,
    keep_multi_map: bool,
    bwa_min_score: Optional[int],
    save_unaligned: bool,
    macs_fdr: Optional[float],
    macs_pvalue: Optional[float],
    save_macs_pileup: bool,
    skip_peak_qc: bool,
    skip_peak_annotation: bool,
    skip_consensus_peaks: bool,
    skip_deseq2_qc: bool,
    skip_fastqc: bool,
    skip_picard_metrics: bool,
    skip_plot_profile: bool,
    skip_plot_fingerprint: bool,
    skip_igv: bool,
    skip_multiqc: bool,
    skip_qc: bool,
    multiqc_methods_description: Optional[LatchFile],
    fragment_size: int = 200,
    narrow_peak: bool = False,
    keep_mito: bool = False,
    min_trimmed_reads: int = 10000,
    aligner: Aligner = Aligner.bwa,
    skip_merge_replicates: bool = False,
    broad_cutoff: float = 0.1,
    min_reps_consensus: int = 1,
    deseq2_vst: bool = True,
    skip_preseq: bool = True,
    skip_ataqv: bool = False,
    save_align_intermeds: bool = True,
    outdir: LatchOutputDir = LatchOutputDir("latch:///Atacseq"),
) -> LatchOutputDir:
    """
    nfcore/atacseq is a bioinformatics analysis pipeline used for ATAC-seq data. LatchBio has made changes to improve performance and created additional downstream outputs.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.2634132-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.2634132)

    **nfcore/atacseq** is a bioinformatics analysis pipeline used for ATAC-seq data.

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    ## Pipeline summary

    1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
    3. Choice of multiple aligners
        1. [`BWA`](https://sourceforge.net/projects/bio-bwa/files/)
        2. [`Chromap`](https://github.com/haowenz/chromap). **For paired-end reads only working until mapping steps, see [here](https://github.com/nf-core/chipseq/issues/291)**
        3. [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
        4. [`STAR`](https://github.com/alexdobin/STAR)
    4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
        1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
        2. Filtering to remove:
            - reads mapping to mitochondrial DNA ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
            - reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
            - reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
            - reads that are not marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
            - reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
            - reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
            - reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
            - reads that are soft-clipped ([`BAMTools`](https://github.com/pezmaster31/bamtools))
            - reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); _paired-end only_)
            - reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
            - reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
            - reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
        3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
        4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
        5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
        6. Calculate genome-wide enrichment (optionally relative to control) ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
        7. Call broad/narrow peaks ([`MACS2`](https://github.com/macs3-project/MACS))
        8. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
        9. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
        10. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
        11. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
        12. Generate ATAC-seq specific QC html report ([`ataqv`](https://github.com/ParkerLab/ataqv))
    6. Merge filtered alignments across replicates ([`picard`](https://broadinstitute.github.io/picard/))
        1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
        2. Remove duplicate reads ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        3. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
        4. Call broad/narrow peaks ([`MACS2`](https://github.com/macs3-project/MACS))
        5. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
        6. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
        7. Count reads in consensus peaks relative to merged library-level alignments ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
        8. Differential accessibility analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    7. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
    8. Present QC for raw read, alignment, peak-calling and differential accessibility results ([`ataqv`](https://github.com/ParkerLab/ataqv), [`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

    ## Usage

    To run on your data, prepare a samplesheet with your input data. Please follow the [documentation on samplesheets](https://nf-co.re/atacseq/usage#samplesheet-input) for more details. An example samplesheet for running the pipeline looks as follows:

    sample,fastq_1,fastq_2,replicate
    CONTROL,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,1
    CONTROL,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,2
    CONTROL,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,3

    See [usage docs](https://nf-co.re/atacseq/usage) for all of the available options when running the pipeline.

    For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/atacseq/usage) and the [parameter documentation](https://nf-co.re/atacseq/parameters).

    ## Pipeline output

    To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/atacseq/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/atacseq/output).

    ## Credits

    The pipeline was originally written by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [Seqera Labs, Spain](https://seqera.io/) and converted to Nextflow DSL2 by Björn Langer ([@bjlang](https://github.com/bjlang)) and Jose Espinosa-Carrasco ([@JoseEspinosa](https://github.com/JoseEspinosa)) from [The Comparative Bioinformatics Group](https://www.crg.eu/en/cedric_notredame) at [The Centre for Genomic Regulation, Spain](https://www.crg.eu/) under the umbrella of the [BovReg project](https://www.bovreg.eu/).
    Many thanks to others who have helped out and contributed along the way too, including (but not limited to): [@ewels](https://github.com/ewels), [@apeltzer](https://github.com/apeltzer), [@crickbabs](https://github.com/crickbabs), [drewjbeh](https://github.com/drewjbeh), [@houghtos](https://github.com/houghtos), [@jinmingda](https://github.com/jinmingda), [@ktrns](https://github.com/ktrns), [@MaxUlysse](https://github.com/MaxUlysse), [@mashehu](https://github.com/mashehu), [@micans](https://github.com/micans), [@pditommaso](https://github.com/pditommaso) and [@sven1103](https://github.com/sven1103).

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).
    For further information or help, don't hesitate to get in touch on the [Slack `#atacseq` channel](https://nfcore.slack.com/channels/atacseq) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/atacseq for your analysis, please cite it using the following doi: [10.5281/zenodo.2634132](https://doi.org/10.5281/zenodo.2634132)
    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
    """

    pvc_name: str = initialize(run_name=run_name)
    NF_Run_Flag = nextflow_runtime(
        pvc_name=pvc_name,
        input=input,
        run_name=run_name,
        fragment_size=fragment_size,
        seq_center=seq_center,
        read_length=read_length,
        with_control=with_control,
        outdir=outdir,
        email=email,
        multiqc_title=multiqc_title,
        genome_source=genome_source,
        genome=genome,
        latch_genome=latch_genome,
        fasta=fasta,
        gtf=gtf,
        gff=gff,
        bwa_index=bwa_index,
        bowtie2_index=bowtie2_index,
        chromap_index=chromap_index,
        star_index=star_index,
        gene_bed=gene_bed,
        tss_bed=tss_bed,
        macs_gsize=macs_gsize,
        blacklist=blacklist,
        mito_name=mito_name,
        save_reference=save_reference,
        keep_mito=keep_mito,
        ataqv_mito_reference=ataqv_mito_reference,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        trim_nextseq=trim_nextseq,
        min_trimmed_reads=min_trimmed_reads,
        skip_trimming=skip_trimming,
        save_trimmed=save_trimmed,
        aligner=aligner,
        keep_dups=keep_dups,
        keep_multi_map=keep_multi_map,
        bwa_min_score=bwa_min_score,
        skip_merge_replicates=skip_merge_replicates,
        save_align_intermeds=save_align_intermeds,
        save_unaligned=save_unaligned,
        narrow_peak=narrow_peak,
        broad_cutoff=broad_cutoff,
        macs_fdr=macs_fdr,
        macs_pvalue=macs_pvalue,
        min_reps_consensus=min_reps_consensus,
        save_macs_pileup=save_macs_pileup,
        skip_peak_qc=skip_peak_qc,
        skip_peak_annotation=skip_peak_annotation,
        skip_consensus_peaks=skip_consensus_peaks,
        deseq2_vst=deseq2_vst,
        skip_deseq2_qc=skip_deseq2_qc,
        skip_fastqc=skip_fastqc,
        skip_picard_metrics=skip_picard_metrics,
        skip_preseq=skip_preseq,
        skip_plot_profile=skip_plot_profile,
        skip_plot_fingerprint=skip_plot_fingerprint,
        skip_igv=skip_igv,
        skip_multiqc=skip_multiqc,
        skip_qc=skip_qc,
        skip_ataqv=skip_ataqv,
        multiqc_methods_description=multiqc_methods_description,
    )
    # NF_Run_Flag = run_name
    input_obj_list, outdir_r_plots = pp.Prepare_Inputs_ATACQC(
        run_flag=NF_Run_Flag, outdir=outdir, aligner=aligner, genome=latch_genome
    )

    input_obj_cov_list, outdir_cov_plots = pp.Prepare_Inputs_Coverages(
        run_flag=NF_Run_Flag, input_dir=outdir, aligner=aligner
    )

    map_tasks_cov = map_task(pp.Compress_Coverages_Sample)(IM=input_obj_cov_list)

    map_task_op = map_task(pp.Run_Rscript)(map_input=input_obj_list)

    Plotting_Data = pp.Calculate_Plotting_Data(
        rplots=map_task_op, cov_plots=map_tasks_cov
    )

    registry_table = pp.UpdateRegistry(d=Plotting_Data, run_name=run_name)

    return outdir


LaunchPlan(
    nf_nf_core_atacseq,
    "Test Data",
    {
        "input": [
            SampleSheet(
                sample="Sample_6",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_6/Sample_6.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_6/Sample_6.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_6",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_6/Sample_6.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_6/Sample_6.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
            SampleSheet(
                sample="Sample_1",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_1/Sample_1.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_1/Sample_1.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_1",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_1/Sample_1.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_1/Sample_1.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
            SampleSheet(
                sample="Sample_2",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_2/Sample_2.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_2/Sample_2.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_2",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_2/Sample_2.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_2/Sample_2.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
            SampleSheet(
                sample="Sample_3",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_3/Sample_3.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_3/Sample_3.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_3",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_3/Sample_3.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_3/Sample_3.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
            SampleSheet(
                sample="Sample_4",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_4/Sample_4.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_4/Sample_4.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_4",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_4/Sample_4.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_4/Sample_4.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
            SampleSheet(
                sample="Sample_5",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_5/Sample_5.Rep_1.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_5/Sample_5.Rep_1.R2.fastq.gz"
                ),
                replicate=1,
            ),
            SampleSheet(
                sample="Sample_5",
                fastq_1=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_5/Sample_5.Rep_2.R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/test-data/35929/Test_Dataset_Verified/Sample_5/Sample_5.Rep_2.R2.fastq.gz"
                ),
                replicate=2,
            ),
        ],
        "genome_source": "input_ref",
        "run_name": "Test-1",
        "read_length": 50,
        "aligner": Aligner.bowtie2,
        "fasta": LatchFile("s3://latch-public/test-data/35929/hg19/hg19.fa"),
        "gtf": LatchFile("s3://latch-public/test-data/35929/hg19/hg19.refGene.gtf"),
        "outdir": LatchOutputDir("latch:///ATAC_Seq_Test"),
    },
)
