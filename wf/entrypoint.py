import os
import shutil
import subprocess
import sys
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import requests
import typing_extensions
from dataclasses_json import dataclass_json
from flytekit.core.annotation import FlyteAnnotation
from latch import map_task, medium_task, small_task
from latch.executions import report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

# import latch_metadata.Compress_Coverage
# from latch_metadata.Compress_Coverage import Compute_Coverage_Across_Samples, Return_Chromosome_Level_Coverages
meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata
import latch_metadata.Compress_Coverage


@dataclass_json
@dataclass
class InputMap:
    run_flag: str
    sample: str
    bamfile: LatchFile
    baifile: LatchFile
    outdir: LatchDir


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    replicate: int


class Reference(Enum):
    hg19 = "GRCh37 (Homo Sapiens hg19)"
    hg38 = "GRCh38 (Homo Sapiens hg38)"


class Aligner(Enum):
    bwa = "bwa"
    bowtie2 = "bowtie2"
    star = "star"
    chromap = "chromap"


class ReadLength(Enum):
    r_50 = "50"
    r_75 = "75"
    r_100 = "100"
    r_150 = "150"
    r_200 = "200"


def latch_listdir(data_directory):
    files = []
    for x in LPath(data_directory).iterdir():
        files.append(x.path)
    return files


def get_flag_defaults(
    name: str, val: typing.Any, default_val: typing.Optional[typing.Any]
):
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)


@custom_task(cpu=4, memory=48, storage_gib=50)
def Calculate_Plotting_Data(
    input_dir: LatchDir, aligner: typing.Optional[Aligner], output_dir: LatchDir
) -> LatchDir:
    data_dir = os.path.join(
        input_dir.local_path, aligner.name, "merged_library/bigwig/"
    )
    local_dir = "cov_parquet/"
    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)
    local_dir = latch_metadata.Compress_Coverage.Compute_Coverage_Across_Samples(
        data_dir, local_dir
    )
    print(local_dir)
    print(os.listdir(local_dir))
    print(os.path.join(output_dir.remote_directory, local_dir))

    return LatchDir(local_dir, os.path.join(output_dir.remote_directory, local_dir))


@custom_task(cpu=4, memory=48, storage_gib=50)
def Run_Rscript(map_input: InputMap) -> LatchDir:
    run_flag = map_input.run_flag
    sample = map_input.sample
    bamfile = map_input.bamfile
    baifile = map_input.baifile
    outdir = map_input.outdir

    local_dir = sample + "/"
    shutil.copy(bamfile.local_path, str(Path().resolve()))
    shutil.copy(baifile.local_path, str(Path().resolve()))

    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)
    bamfile.download()
    baifile.download()

    cmd_RunRscript = [
        "mamba",
        "run",
        "-n",
        "atacseqqc",
        "Rscript",
        "/root/latch_metadata/ATACSeqQC_Plots.R",
        Path(bamfile.local_path).name,
        local_dir,
    ]
    print(" ".join(cmd_RunRscript))
    subprocess.run(" ".join(cmd_RunRscript), shell=True, check=True)
    return LatchDir(local_dir, os.path.join(outdir.remote_directory, sample))


@small_task
def Prepare_Inputs(
    run_flag: str, outdir: LatchDir, aligner: typing.Optional[Aligner] = Aligner.bwa
) -> (typing.List[InputMap], LatchDir):
    input_path = os.path.join(outdir.remote_directory, aligner.name, "merged_library")
    files = latch_listdir(input_path)
    if not os.path.isdir("R_Plots/"):
        os.mkdir("R_Plots/")
    outdir_r_plots = LatchDir(
        "R_Plots", os.path.join(outdir.remote_directory, "R_Plots")
    )
    object_list = []
    for f in files:
        if f.endswith(".mLb.mkD.sorted.bam"):
            bamfile = f
            baifile = f + ".bai"
            assert baifile in files, "Missing bai file for " + baifile
            sample = f.split("/")[-1].replace(".mLb.mkD.sorted.bam", "")
            o = InputMap(run_flag, sample, bamfile, baifile, outdir_r_plots)
            object_list.append(o)
    return object_list, outdir_r_plots


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


input_construct_samplesheet = metadata._nextflow_metadata.parameters[
    "input"
].samplesheet_constructor


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    pvc_name: str,
    input: typing.List[SampleSheet],
    seq_center: typing.Optional[str],
    read_length: typing.Optional[ReadLength],
    with_control: typing.Optional[bool],
    outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
    email: typing.Optional[str],
    multiqc_title: typing.Optional[str],
    genome_source: str,
    genome: typing.Optional[Reference],
    fasta: typing.Optional[LatchFile],
    gtf: typing.Optional[LatchFile],
    gff: typing.Optional[LatchFile],
    bwa_index: typing.Optional[str],
    bowtie2_index: typing.Optional[LatchDir],
    chromap_index: typing.Optional[str],
    star_index: typing.Optional[str],
    gene_bed: typing.Optional[LatchFile],
    tss_bed: typing.Optional[LatchFile],
    macs_gsize: typing.Optional[float],
    blacklist: typing.Optional[str],
    mito_name: typing.Optional[str],
    save_reference: typing.Optional[bool],
    ataqv_mito_reference: typing.Optional[str],
    clip_r1: typing.Optional[int],
    clip_r2: typing.Optional[int],
    three_prime_clip_r1: typing.Optional[int],
    three_prime_clip_r2: typing.Optional[int],
    trim_nextseq: typing.Optional[int],
    skip_trimming: typing.Optional[bool],
    save_trimmed: typing.Optional[bool],
    keep_dups: typing.Optional[bool],
    keep_multi_map: typing.Optional[bool],
    bwa_min_score: typing.Optional[int],
    save_align_intermeds: typing.Optional[bool],
    save_unaligned: typing.Optional[bool],
    narrow_peak: typing.Optional[bool],
    macs_fdr: typing.Optional[float],
    macs_pvalue: typing.Optional[float],
    save_macs_pileup: typing.Optional[bool],
    skip_peak_qc: typing.Optional[bool],
    skip_peak_annotation: typing.Optional[bool],
    skip_consensus_peaks: typing.Optional[bool],
    skip_deseq2_qc: typing.Optional[bool],
    skip_fastqc: typing.Optional[bool],
    skip_picard_metrics: typing.Optional[bool],
    skip_plot_profile: typing.Optional[bool],
    skip_plot_fingerprint: typing.Optional[bool],
    skip_igv: typing.Optional[bool],
    skip_multiqc: typing.Optional[bool],
    skip_qc: typing.Optional[bool],
    multiqc_methods_description: typing.Optional[LatchFile],
    fragment_size: typing.Optional[int],
    keep_mito: typing.Optional[bool],
    min_trimmed_reads: typing.Optional[int],
    aligner: typing.Optional[Aligner],
    skip_merge_replicates: typing.Optional[bool],
    broad_cutoff: typing.Optional[float],
    min_reps_consensus: typing.Optional[int],
    deseq2_vst: typing.Optional[bool],
    skip_preseq: typing.Optional[bool],
    skip_ataqv: typing.Optional[bool],
) -> str:
    shared_dir = Path("/nf-workdir")

    input_samplesheet = input_construct_samplesheet(input)

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    profile_list = ["docker"]
    if False:
        profile_list.extend([p.value for p in execution_profiles])

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ",".join(profile_list)

    if genome_source == "latch_genome_source":
        fasta = os.path.join(
            "s3://latch-public/test-data/35929/", genome.name, genome.name + ".fa"
        )
        gtf = os.path.join(
            "s3://latch-public/test-data/35929/",
            genome.name,
            genome.name + ".refGene.gtf",
        )
        if aligner.name == "bowtie2":
            bowtie2_index = os.path.join(
                "s3://latch-public/test-data/35929/", genome.name, "bowtie2"
            )
        elif aligner.name == "bwa":
            bwa_index = os.path.join(
                "s3://latch-public/test-data/35929/", genome.name, "bwa"
            )
    print(aligner.name, genome, fasta, bowtie2_index)

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        profiles,
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input_samplesheet),
        *get_flag_defaults("fragment_size", fragment_size, 200),
        *get_flag_defaults("seq_center", seq_center, None),
        *get_flag_defaults("read_length", read_length, None),
        *get_flag_defaults("with_control", with_control, None),
        *get_flag("outdir", outdir),
        *get_flag_defaults("email", email, None),
        *get_flag_defaults("multiqc_title", multiqc_title, None),
        # *get_flag_defaults("genome", genome, None),
        *get_flag_defaults("fasta", fasta, None),
        *get_flag_defaults("gtf", gtf, None),
        *get_flag_defaults("gff", gff, None),
        *get_flag_defaults("bwa_index", bwa_index, None),
        *get_flag_defaults("bowtie2_index", bowtie2_index, None),
        *get_flag_defaults("chromap_index", chromap_index, None),
        *get_flag_defaults("star_index", star_index, None),
        *get_flag_defaults("gene_bed", gene_bed, None),
        *get_flag_defaults("tss_bed", tss_bed, None),
        *get_flag_defaults("macs_gsize", macs_gsize, None),
        *get_flag_defaults("blacklist", blacklist, None),
        *get_flag_defaults("mito_name", mito_name, None),
        *get_flag_defaults("save_reference", save_reference, None),
        *get_flag_defaults("keep_mito", keep_mito, None),
        *get_flag_defaults("ataqv_mito_reference", ataqv_mito_reference, None),
        *get_flag_defaults("clip_r1", clip_r1, None),
        *get_flag_defaults("clip_r2", clip_r2, None),
        *get_flag_defaults("three_prime_clip_r1", three_prime_clip_r1, None),
        *get_flag_defaults("three_prime_clip_r2", three_prime_clip_r2, None),
        *get_flag_defaults("trim_nextseq", trim_nextseq, None),
        *get_flag_defaults("min_trimmed_reads", min_trimmed_reads, 10000),
        *get_flag_defaults("skip_trimming", skip_trimming, None),
        *get_flag_defaults("save_trimmed", save_trimmed, None),
        *get_flag_defaults("aligner", aligner, Aligner.bwa),
        *get_flag_defaults("keep_dups", keep_dups, None),
        *get_flag_defaults("keep_multi_map", keep_multi_map, None),
        *get_flag_defaults("bwa_min_score", bwa_min_score, None),
        *get_flag_defaults("skip_merge_replicates", skip_merge_replicates, None),
        *get_flag_defaults("save_align_intermeds", save_align_intermeds, None),
        *get_flag_defaults("save_unaligned", save_unaligned, None),
        *get_flag_defaults("narrow_peak", narrow_peak, None),
        *get_flag_defaults("broad_cutoff", broad_cutoff, 0.1),
        *get_flag_defaults("macs_fdr", macs_fdr, None),
        *get_flag_defaults("macs_pvalue", macs_pvalue, None),
        *get_flag_defaults("min_reps_consensus", min_reps_consensus, 1),
        *get_flag_defaults("save_macs_pileup", save_macs_pileup, None),
        *get_flag_defaults("skip_peak_qc", skip_peak_qc, None),
        *get_flag_defaults("skip_peak_annotation", skip_peak_annotation, None),
        *get_flag_defaults("skip_consensus_peaks", skip_consensus_peaks, None),
        *get_flag_defaults("deseq2_vst", deseq2_vst, True),
        *get_flag_defaults("skip_deseq2_qc", skip_deseq2_qc, None),
        *get_flag_defaults("skip_fastqc", skip_fastqc, None),
        *get_flag_defaults("skip_picard_metrics", skip_picard_metrics, None),
        *get_flag_defaults("skip_preseq", skip_preseq, True),
        *get_flag_defaults("skip_plot_profile", skip_plot_profile, None),
        *get_flag_defaults("skip_plot_fingerprint", skip_plot_fingerprint, None),
        *get_flag_defaults("skip_igv", skip_igv, None),
        *get_flag_defaults("skip_multiqc", skip_multiqc, None),
        *get_flag_defaults("skip_qc", skip_qc, None),
        *get_flag_defaults("skip_ataqv", skip_ataqv, None),
        *get_flag_defaults(
            "multiqc_methods_description", multiqc_methods_description, None
        ),
    ]

    print("Launching Nextflow Runtime")
    print(cmd)
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///your_log_dir/nf_nf_core_atacseq", name, "nextflow.log"
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)
    return "Done"


@workflow(metadata._nextflow_metadata)
def nf_nf_core_atacseq(
    input: typing.List[SampleSheet],
    seq_center: typing.Optional[str],
    read_length: typing.Optional[ReadLength],
    with_control: typing.Optional[bool],
    outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
    email: typing.Optional[str],
    multiqc_title: typing.Optional[str],
    genome_source: str,
    genome: typing.Optional[Reference],
    fasta: typing.Optional[LatchFile],
    gtf: typing.Optional[LatchFile],
    gff: typing.Optional[LatchFile],
    bwa_index: typing.Optional[str],
    bowtie2_index: typing.Optional[LatchDir],
    chromap_index: typing.Optional[str],
    star_index: typing.Optional[str],
    gene_bed: typing.Optional[LatchFile],
    tss_bed: typing.Optional[LatchFile],
    macs_gsize: typing.Optional[float],
    blacklist: typing.Optional[str],
    mito_name: typing.Optional[str],
    save_reference: typing.Optional[bool],
    ataqv_mito_reference: typing.Optional[str],
    clip_r1: typing.Optional[int],
    clip_r2: typing.Optional[int],
    three_prime_clip_r1: typing.Optional[int],
    three_prime_clip_r2: typing.Optional[int],
    trim_nextseq: typing.Optional[int],
    skip_trimming: typing.Optional[bool],
    save_trimmed: typing.Optional[bool],
    keep_dups: typing.Optional[bool],
    keep_multi_map: typing.Optional[bool],
    bwa_min_score: typing.Optional[int],
    save_align_intermeds: typing.Optional[bool],
    save_unaligned: typing.Optional[bool],
    narrow_peak: typing.Optional[bool],
    macs_fdr: typing.Optional[float],
    macs_pvalue: typing.Optional[float],
    save_macs_pileup: typing.Optional[bool],
    skip_peak_qc: typing.Optional[bool],
    skip_peak_annotation: typing.Optional[bool],
    skip_consensus_peaks: typing.Optional[bool],
    skip_deseq2_qc: typing.Optional[bool],
    skip_fastqc: typing.Optional[bool],
    skip_picard_metrics: typing.Optional[bool],
    skip_plot_profile: typing.Optional[bool],
    skip_plot_fingerprint: typing.Optional[bool],
    skip_igv: typing.Optional[bool],
    skip_multiqc: typing.Optional[bool],
    skip_qc: typing.Optional[bool],
    multiqc_methods_description: typing.Optional[LatchFile],
    fragment_size: typing.Optional[int] = 200,
    keep_mito: typing.Optional[bool] = False,
    min_trimmed_reads: typing.Optional[int] = 10000,
    aligner: typing.Optional[Aligner] = Aligner.bwa,
    skip_merge_replicates: typing.Optional[bool] = False,
    broad_cutoff: typing.Optional[float] = 0.1,
    min_reps_consensus: typing.Optional[int] = 1,
    deseq2_vst: typing.Optional[bool] = True,
    skip_preseq: typing.Optional[bool] = True,
    skip_ataqv: typing.Optional[bool] = False,
) -> LatchDir:
    """
    nf-core/atacseq

    Sample Description
    """

    pvc_name: str = initialize()
    NF_Run_Flag = nextflow_runtime(
        pvc_name=pvc_name,
        input=input,
        fragment_size=fragment_size,
        seq_center=seq_center,
        read_length=read_length,
        with_control=with_control,
        outdir=outdir,
        email=email,
        multiqc_title=multiqc_title,
        genome_source=genome_source,
        genome=genome,
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
    input_obj_list, outdir_r_plots = Prepare_Inputs(
        run_flag=NF_Run_Flag, outdir=outdir, aligner=aligner
    )

    map_task_op = map_task(Run_Rscript)(map_input=input_obj_list)

    outdir_cov = Calculate_Plotting_Data(
        input_dir=outdir, aligner=aligner, output_dir=outdir_r_plots
    )
    # print(outdir_cov)

    return outdir
