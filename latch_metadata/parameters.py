
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

generated_parameters = {
    'input': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title='Input/output options',
        description='Path to comma-separated file containing information about the samples in the experiment.',
    ),
    'fragment_size': NextflowParameter(
        type=typing.Optional[int],
        default=200,
        section_title=None,
        description='Estimated fragment size used to extend single-end reads.',
    ),
    'seq_center': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Sequencing center information to be added to read group of BAM files.',
    ),
    'read_length': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description="Read length used to calculate or retrieve pre-computed MACS2 genome size for peak calling if `--macs_gsize` isn't provided.",
    ),
    'with_control': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Use controls.',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'email': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Email address for completion summary.',
    ),
    'multiqc_title': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='MultiQC report title. Printed as page header, used for filename if not otherwise specified.',
    ),
    'genome': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title='Reference genome options',
        description='Name of iGenomes reference.',
    ),
    'fasta': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to FASTA genome file.',
    ),
    'gtf': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to GTF annotation file.',
    ),
    'gff': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to GFF3 annotation file.',
    ),
    'bwa_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to directory or tar.gz archive for pre-built BWA index.',
    ),
    'bowtie2_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to directory or tar.gz archive for pre-built Bowtie2 index.',
    ),
    'chromap_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to directory or tar.gz archive for pre-built Chromap index.',
    ),
    'star_index': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to directory or tar.gz archive for pre-built STAR index.',
    ),
    'gene_bed': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to BED file containing gene intervals. This will be created from the GTF file if not specified.',
    ),
    'tss_bed': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to BED file containing transcription start sites. This will be created from the gene BED file if not specified.',
    ),
    'macs_gsize': NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        description='Effective genome size parameter required by MACS2.',
    ),
    'blacklist': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to blacklist regions in BED format, used for filtering alignments.',
    ),
    'mito_name': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Name of Mitochondrial chomosome in reference assembly e.g. chrM.',
    ),
    'save_reference': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='If generated by the pipeline save the aligner index (e.g. BWA) in the results directory.',
    ),
    'keep_mito': NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        description='Reads mapping to mitochondrial contig are not filtered from alignments.',
    ),
    'ataqv_mito_reference': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Sets the value of the ataqv --mitochondrial-reference-name argument',
    ),
    'clip_r1': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title='Adapter trimming options',
        description="Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).",
    ),
    'clip_r2': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).",
    ),
    'three_prime_clip_r1': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.",
    ),
    'three_prime_clip_r2': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.",
    ),
    'trim_nextseq': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails.',
    ),
    'min_trimmed_reads': NextflowParameter(
        type=typing.Optional[int],
        default=10000,
        section_title=None,
        description='Minimum number of trimmed reads below which samples are removed from further processing. Some downstream steps in the pipeline will fail if this threshold is too low.',
    ),
    'skip_trimming': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip the adapter trimming step.',
    ),
    'save_trimmed': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save the trimmed FastQ files in the results directory.',
    ),
    'aligner': NextflowParameter(
        type=typing.Optional[str],
        default='bwa',
        section_title='Alignment options',
        description="Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2', 'chromap' and 'star'.",
    ),
    'keep_dups': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Duplicate reads are not filtered from alignments.',
    ),
    'keep_multi_map': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Reads mapping to multiple locations are not filtered from alignments.',
    ),
    'bwa_min_score': NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        description='Don’t output BWA MEM alignments with score lower than this parameter.',
    ),
    'skip_merge_replicates': NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        description='Do not perform alignment merging and downstream analysis by merging replicates i.e. only do this by merging resequenced libraries.',
    ),
    'save_align_intermeds': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save the intermediate BAM files from the alignment step.',
    ),
    'save_unaligned': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory.',
    ),
    'narrow_peak': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Peak calling options',
        description='Run MACS2 in narrowPeak mode.',
    ),
    'broad_cutoff': NextflowParameter(
        type=typing.Optional[float],
        default=0.1,
        section_title=None,
        description='Specifies broad cutoff value for MACS2. Only used when --narrow_peak isnt specified.',
    ),
    'macs_fdr': NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        description='Minimum FDR (q-value) cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive.',
    ),
    'macs_pvalue': NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        description='p-value cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive. If --macs_pvalue cutoff is set, q-value will not be calculated and reported as -1 in the final .xls file.',
    ),
    'min_reps_consensus': NextflowParameter(
        type=typing.Optional[int],
        default=1,
        section_title=None,
        description='Number of biological replicates required from a given condition for a peak to contribute to a consensus peak.',
    ),
    'save_macs_pileup': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Instruct MACS2 to create bedGraph files normalised to signal per million reads.',
    ),
    'skip_peak_qc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip MACS2 peak QC plot generation.',
    ),
    'skip_peak_annotation': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip annotation of MACS2 and consensus peaks with HOMER.',
    ),
    'skip_consensus_peaks': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip consensus peak generation, annotation and counting.',
    ),
    'deseq2_vst': NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title='Differential analysis options',
        description='Use vst transformation instead of rlog with DESeq2.',
    ),
    'skip_deseq2_qc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip DESeq2 PCA and heatmap plotting.',
    ),
    'skip_fastqc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Process skipping options',
        description='Skip FastQC.',
    ),
    'skip_picard_metrics': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip Picard CollectMultipleMetrics.',
    ),
    'skip_preseq': NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title=None,
        description='Skip Preseq.',
    ),
    'skip_plot_profile': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip deepTools plotProfile.',
    ),
    'skip_plot_fingerprint': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip deepTools plotFingerprint.',
    ),
    'skip_igv': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip IGV.',
    ),
    'skip_multiqc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip MultiQC.',
    ),
    'skip_qc': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Skip all QC steps except for MultiQC.',
    ),
    'skip_ataqv': NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        description='Skip Ataqv.',
    ),
    'multiqc_methods_description': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title='Generic options',
        description='Custom MultiQC yaml file containing HTML including a methods description.',
    ),
}

