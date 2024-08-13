import csv
import typing
from dataclasses import dataclass
from enum import Enum

import typing_extensions
from flytekit.core.annotation import FlyteAnnotation
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import NextflowParameter


# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters
class Aligner(Enum):
    bwa = "bwa"
    bowtie2 = "bowtie2"
    star = "star"
    chromap = "chromap"


class Reference(Enum):
    hg19 = "GRCh37 (Homo Sapiens hg19)"
    hg38 = "GRCh38 (Homo Sapiens hg38)"


class ReadLength(Enum):
    r_50 = "50"
    r_75 = "75"
    r_100 = "100"
    r_150 = "150"
    r_200 = "200"


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    replicate: int


generated_parameters = {
    "input": NextflowParameter(
        type=typing.List[SampleSheet],
        display_name="Samplesheet",
        description="Information about the samples in the experiment",
        section_title=None,
        samplesheet_type="csv",
        samplesheet=True,
    ),
    "fragment_size": NextflowParameter(
        type=typing.Optional[int],
        default=200,
        display_name="Fragment Size",
        section_title=None,
        description="Estimated fragment size used to extend single-end reads.",
    ),
    "seq_center": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        display_name="Sequencing Center",
        section_title=None,
        description="Sequencing center information to be added to read group of BAM files.",
    ),
    "read_length": NextflowParameter(
        type=typing.Optional[ReadLength],
        default=None,
        section_title=None,
        display_name="Read Length",
        description="Read length used to calculate or retrieve pre-computed MACS2 genome size for peak calling if `--macs_gsize` isn't provided.",
    ),
    "with_control": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Use Controls",
        description="Use controls.",
    ),
    "outdir": NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
        default=None,
        section_title=None,
        display_name="Output Directory",
        description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
    ),
    "email": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        display_name="email",
        section_title=None,
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Multi QC Title",
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "genome_source": NextflowParameter(
        type=str,
        display_name="Reference Genome",
        description="Choose Reference Genome",
    ),
    "genome": NextflowParameter(
        type=typing.Optional[Reference],
        default=None,
        section_title="Reference genome options",
        display_name="Genome",
        description="Name of the Latch Verfied Reference Genome",
    ),
    "fasta": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Fasta File",
        description="Path to FASTA genome file.",
    ),
    "gtf": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="GTF",
        description="Path to GTF annotation file.",
    ),
    "gff": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="GFF3",
        description="Path to GFF3 annotation file.",
    ),
    "bwa_index": NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        display_name="BWA Index",
        description="Path to directory or tar.gz archive for pre-built BWA index.",
    ),
    "bowtie2_index": NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        display_name="Bowtie2 Index",
        description="Path to directory or tar.gz archive for pre-built Bowtie2 index.",
    ),
    "chromap_index": NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        display_name="Chromap Index",
        description="Path to directory or tar.gz archive for pre-built Chromap index.",
    ),
    "star_index": NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        display_name="Star Index",
        description="Path to directory or tar.gz archive for pre-built STAR index.",
    ),
    "gene_bed": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="Gene Bed File",
        description="Path to BED file containing gene intervals. This will be created from the GTF file if not specified.",
    ),
    "tss_bed": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        display_name="TSS Bed File",
        description="Path to BED file containing transcription start sites. This will be created from the gene BED file if not specified.",
    ),
    "macs_gsize": NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        display_name="Genome Size",
        description="Effective genome size parameter required by MACS2. Either the read length or the macs_gsize must be provided.",
    ),
    "blacklist": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Blacklisted Regions",
        description="Path to blacklist regions in BED format, used for filtering alignments.",
    ),
    "mito_name": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="Name of the Mito. Chromsome",
        description="Name of Mitochondrial chomosome in reference assembly e.g. chrM.",
    ),
    "save_reference": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Save Reference",
        description="If generated by the pipeline save the aligner index (e.g. BWA) in the results directory.",
    ),
    "keep_mito": NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        display_name="Keep Mitrochondial Region",
        description="Reads mapping to mitochondrial contig are not filtered from alignments.",
    ),
    "ataqv_mito_reference": NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        display_name="--mitochondrial-reference-name",
        description="Sets the value of the ataqv --mitochondrial-reference-name argument",
    ),
    "clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title="Adapter trimming options",
        display_name="Length of the 5' end to be clipped of foward read",
        description="Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).",
    ),
    "clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        display_name="Length of the 5' end to be clipped of reverse read",
        description="Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        display_name="Length of the 3' end to be clipped of the forward read",
        description="Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        display_name="Length of the 3' end to be clipped of the reverse read",
        description="Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.",
    ),
    "trim_nextseq": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        display_name="Apply nextseq filter",
        description="Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails.",
    ),
    "min_trimmed_reads": NextflowParameter(
        type=typing.Optional[int],
        default=10000,
        section_title=None,
        display_name="Minimum number of trimmed reads",
        description="Minimum number of trimmed reads below which samples are removed from further processing. Some downstream steps in the pipeline will fail if this threshold is too low.",
    ),
    "skip_trimming": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip adapter trimming",
        description="Skip the adapter trimming step.",
    ),
    "save_trimmed": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Save trimmed reads",
        description="Save the trimmed FastQ files in the results directory.",
    ),
    "aligner": NextflowParameter(
        type=typing.Optional[Aligner],
        default=Aligner.bwa,
        section_title="Alignment options",
        display_name="Alignment options",
        description="Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2', 'chromap' and 'star'.",
    ),
    "keep_dups": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="keep duplicates",
        description="Duplicate reads are not filtered from alignments.",
    ),
    "keep_multi_map": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Keep Multimapped Reads",
        description="Reads mapping to multiple locations are not filtered from alignments.",
    ),
    "bwa_min_score": NextflowParameter(
        type=typing.Optional[int],
        default=None,
        section_title=None,
        display_name="Min. Score for BWA Mem",
        description="Don’t output BWA MEM alignments with score lower than this parameter.",
    ),
    "skip_merge_replicates": NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        display_name="Skip Merging Replicates",
        description="Do not perform alignment merging and downstream analysis by merging replicates i.e. only do this by merging resequenced libraries.",
    ),
    "save_align_intermeds": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Save Intermediate Alignments",
        description="Save the intermediate BAM files from the alignment step.",
    ),
    "save_unaligned": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Save Unaligned Reads",
        description="Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory.",
    ),
    "narrow_peak": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title="Peak calling options",
        display_name="Narro Peaks",
        description="Run MACS2 in narrowPeak mode.",
    ),
    "broad_cutoff": NextflowParameter(
        type=typing.Optional[float],
        default=0.1,
        section_title=None,
        display_name="Cutoff for Broad Peaks",
        description="Specifies broad cutoff value for MACS2. Only used when --narrow_peak isnt specified.",
    ),
    "macs_fdr": NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        display_name="MACS2 FDR",
        description="Minimum FDR (q-value) cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive.",
    ),
    "macs_pvalue": NextflowParameter(
        type=typing.Optional[float],
        default=None,
        section_title=None,
        display_name="MACS2 P-Value",
        description="p-value cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive. If --macs_pvalue cutoff is set, q-value will not be calculated and reported as -1 in the final .xls file.",
    ),
    "min_reps_consensus": NextflowParameter(
        type=typing.Optional[int],
        default=1,
        section_title=None,
        display_name="Number of Minum Replicates",
        description="Number of biological replicates required from a given condition for a peak to contribute to a consensus peak.",
    ),
    "save_macs_pileup": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Save MACS2 Pileup",
        description="Instruct MACS2 to create bedGraph files normalised to signal per million reads.",
    ),
    "skip_peak_qc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Peak QC",
        description="Skip MACS2 peak QC plot generation.",
    ),
    "skip_peak_annotation": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Peak Annotation",
        description="Skip annotation of MACS2 and consensus peaks with HOMER.",
    ),
    "skip_consensus_peaks": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Peak Consensus",
        description="Skip consensus peak generation, annotation and counting.",
    ),
    "deseq2_vst": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title="Differential analysis options",
        display_name="VST Transformation with DESeq2",
        description="Use vst transformation instead of rlog with DESeq2.",
    ),
    "skip_deseq2_qc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip DESeq2 QC",
        description="Skip DESeq2 PCA and heatmap plotting.",
    ),
    "skip_fastqc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title="Process skipping options",
        display_name="Skip FastQC",
        description="Skip FastQC.",
    ),
    "skip_picard_metrics": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Picard Metrics",
        description="Skip Picard CollectMultipleMetrics.",
    ),
    "skip_preseq": NextflowParameter(
        type=typing.Optional[bool],
        default=True,
        section_title=None,
        display_name="Skip Preseq",
        description="Skip Preseq.",
    ),
    "skip_plot_profile": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Plot Profile",
        description="Skip deepTools plotProfile.",
    ),
    "skip_plot_fingerprint": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip Plot Fingerprint",
        description="Skip deepTools plotFingerprint.",
    ),
    "skip_igv": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip IGV",
        description="Skip IGV.",
    ),
    "skip_multiqc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip MultiQC",
        description="Skip MultiQC.",
    ),
    "skip_qc": NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        display_name="Skip QC",
        description="Skip all QC steps except for MultiQC.",
    ),
    "skip_ataqv": NextflowParameter(
        type=typing.Optional[bool],
        default=False,
        section_title=None,
        display_name="Skip ATAQV",
        description="Skip Ataqv.",
    ),
    "multiqc_methods_description": NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title="Generic options",
        display_name="MultiQC Description",
        description="Custom MultiQC yaml file containing HTML including a methods description.",
    ),
}
