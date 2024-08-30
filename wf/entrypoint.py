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
from latch.account import Account
from latch.executions import report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.registry.project import Project
from latch.registry.table import Table
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
class InputMap_ATACQC:
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


@dataclass
class Registry_Obj:
    sample: str
    frag_file: LatchFile
    feat_alignment: LatchFile
    sat_curves: LatchFile
    cov_parquet: LatchDir


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
    for x in data_directory.iterdir():
        files.append(x.path)
    return files


def get_flag_defaults(
    name: str, val: typing.Any, default_val: typing.Optional[typing.Any]
):
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)


@dataclass_json
@dataclass
class InputMap_Compress_Coverages:
    data_path: LatchFile
    outPath: LatchDir
    sample: str


@custom_task(cpu=4, memory=48, storage_gib=50)
def Prepare_Inputs_Coverages(
    run_flag: str, input_dir: LatchDir, aligner: typing.Optional[Aligner] = Aligner.bwa
) -> (typing.List[InputMap_Compress_Coverages], LatchDir):
    bigwig_files = LPath(
        os.path.join(
            input_dir.remote_path, run_flag, aligner.name, "merged_library/bigwig/"
        )
    )
    local_dir = "cov_parquet/"
    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)
    outdir_cov_plots = LatchDir(
        local_dir, f"{input_dir.remote_directory}/{run_flag}/{local_dir}"
    )

    input_objects = []
    for b in latch_listdir(bigwig_files):
        if b.endswith("bigWig"):
            sample = Path(b).name
            print(sample)
            o = InputMap_Compress_Coverages(
                LatchFile(b), outdir_cov_plots, sample.replace(".mLb.clN.bigWig", "")
            )
            input_objects.append(o)
    return input_objects, outdir_cov_plots


@custom_task(cpu=4, memory=48, storage_gib=50)
def Merge_pq_files(
    cov_parquet_path: typing.List[LatchDir], outPath: LatchDir
) -> LatchDir:
    chromosomes = {}
    for cov in cov_parquet_path:
        files = latch_listdir(cov)
        for s in files:
            print(s)
            splits = Path(s).name.split(".")
            print(splits[0])
            try:
                chromosomes[splits[0]].append(LatchFile(s))
            except KeyError:
                chromosomes[splits[0]] = [LatchFile(s)]

    out_dir = "merged_parquet_files"
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for c in chromosomes:
        d = {}
        keys = chromosomes[c]
        for k in keys:
            table = latch_metadata.Compress_Coverage.pq.read_table(k.local_path)
            sample = Path(k.local_path).name.replace(".parquet", "")
            print(c, sample)
            d[sample] = table["data"]
        latch_metadata.Compress_Coverage.pq.write_table(
            latch_metadata.Compress_Coverage.pa.table(d),
            os.path.join(out_dir, c + ".parquet"),
        )
    return LatchDir(out_dir, os.path.join(outPath.remote_directory, out_dir))


@custom_task(cpu=4, memory=48, storage_gib=50)
def Compress_Coverages_Sample(IM: InputMap_Compress_Coverages) -> LatchDir:
    bigWig_file = IM.data_path
    sample = IM.sample
    outPath = IM.outPath

    local_dir = sample

    print(bigWig_file, outPath, sample)
    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)

    latch_metadata.Compress_Coverage.Compress_Coverage(
        bigWig_file.local_path, local_dir, sample
    )

    return LatchDir(local_dir, os.path.join(outPath.remote_directory, local_dir))


@small_task
def Calculate_Plotting_Data(
    rplots: typing.List[LatchDir], cov_plots: typing.List[LatchDir]
) -> typing.List[Registry_Obj]:
    d = {}
    for f in rplots:
        files = latch_listdir(f)
        sample = Path(f.remote_path).name
        plots = {}
        for i in files:
            if Path(i).name in [
                "Frag_Sizes.txt",
                "featurealignment_coverage.txt",
                "Saturation_Plots.txt",
            ]:
                plots[Path(i).name] = LatchFile(i)
        d[sample] = plots

    for f in cov_plots:
        sample = Path(f.remote_path).name
        d[sample]["Cov_Plots"] = f

    obj_list = []
    for k in d:
        print(k)
        o = Registry_Obj(
            k,
            d[k]["Frag_Sizes.txt"],
            d[k]["featurealignment_coverage.txt"],
            d[k]["Saturation_Plots.txt"],
            d[k]["Cov_Plots"],
        )
        obj_list.append(o)

    return obj_list


@custom_task(cpu=4, memory=48, storage_gib=50)
def Run_Rscript(map_input: InputMap_ATACQC) -> LatchDir:
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
def UpdateRegistry(d: typing.List[Registry_Obj], run_name: str) -> str:
    target_project_name = "ATAC_Seq_Results"
    target_table_name = f"{run_name}_Results"

    account = Account.current()
    target_project = next(
        (
            project
            for project in account.list_registry_projects()
            if project.get_display_name() == target_project_name
        ),
        None,
    )

    if target_project is None:
        with account.update() as account_updater:
            account_updater.upsert_registry_project(target_project_name)
        target_project = next(
            project
            for project in account.list_registry_projects()
            if project.get_display_name() == target_project_name
        )
        print("Upserted project")

    target_table = next(
        (
            table
            for table in target_project.list_tables()
            if table.get_display_name() == target_table_name
        ),
        None,
    )

    columns = [
        "sample",
        "Fragment_Size_Distribution",
        "Feature_Alignment_Coverage",
        "Saturation_Curves",
        "Cov_Parquet_Files",
    ]
    col_types = [str, LatchFile, LatchFile, LatchFile, LatchDir]
    if target_table == None:
        ### Upsert_Table
        with target_project.update() as project_updater:
            project_updater.upsert_table(target_table_name)
        target_table = next(
            table
            for table in target_project.list_tables()
            if table.get_display_name() == target_table_name
        )
        print("Upserted table")

        with target_table.update() as updater:
            for i in range(0, len(columns)):
                updater.upsert_column(columns[i], type=col_types[i])
        print("Upserted columns")
    ctr = len(target_table.get_dataframe())
    with target_table.update() as updater:
        for v in d:
            updater.upsert_record(
                name=str(ctr),
                sample=v.sample,
                Fragment_Size_Distribution=v.frag_file,
                Feature_Alignment_Coverage=v.feat_alignment,
                Saturation_Curves=v.sat_curves,
                Cov_Parquet_Files=v.cov_parquet,
            )
            ctr += 1

    return target_table.id


@small_task
def Prepare_Inputs_ATACQC(
    run_flag: str, outdir: LatchDir, aligner: typing.Optional[Aligner] = Aligner.bwa
) -> (typing.List[InputMap_ATACQC], LatchDir):
    input_path = LPath(
        os.path.join(outdir.remote_directory, run_flag, aligner.name, "merged_library")
    )
    files = latch_listdir(input_path)
    if not os.path.isdir("R_Plots/"):
        os.mkdir("R_Plots/")
    outdir_r_plots = LatchDir(
        "R_Plots", os.path.join(f"{outdir.remote_directory}/{run_flag}/R_Plots")
    )
    object_list = []
    for f in files:
        print(f)
        if f.endswith(".mLb.mkD.sorted.bam"):
            bamfile = f
            baifile = f + ".bai"
            assert baifile in files, "Missing bai file for " + baifile
            sample = f.split("/")[-1].replace(".mLb.mkD.sorted.bam", "")
            o = InputMap_ATACQC(run_flag, sample, bamfile, baifile, outdir_r_plots)
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
    run_name: str,
    seq_center: typing.Optional[str],
    read_length: typing.Optional[int],
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
    save_align_intermeds: bool,
    save_unaligned: typing.Optional[bool],
    narrow_peak: bool,
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
    fragment_size: int,
    keep_mito: bool,
    min_trimmed_reads: int,
    aligner: Aligner,
    skip_merge_replicates: bool,
    broad_cutoff: float,
    min_reps_consensus: int,
    deseq2_vst: bool,
    skip_preseq: bool,
    skip_ataqv: bool,
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
        *get_flag_defaults(
            "outdir", LatchOutputDir(outdir.remote_path + "/" + run_name), None
        ),
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
    return run_name


@workflow(metadata._nextflow_metadata)
def nf_nf_core_atacseq(
    input: typing.List[SampleSheet],
    run_name: str,
    seq_center: typing.Optional[str],
    read_length: typing.Optional[int],
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
    save_unaligned: typing.Optional[bool],
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
) -> LatchDir:
    """
    Latch Verified ATAC Seq

    # Latch Verified ATAC Seq
    **nfcore/atacseq** is a bioinformatics analysis pipeline used for ATAC-seq data.

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!
    On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/atacseq/results).

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

    > **Note**
    > If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
    > to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
    > with `-profile test` before running the workflow on actual data.

    To run on your data, prepare a tab-separated samplesheet with your input data. Please follow the [documentation on samplesheets](https://nf-co.re/atacseq/usage#samplesheet-input) for more details. An example samplesheet for running the pipeline looks as follows:

    sample,fastq_1,fastq_2,replicate
    CONTROL,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,1
    CONTROL,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,2
    CONTROL,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,3

    See [usage docs](https://nf-co.re/atacseq/usage) for all of the available options when running the pipeline.

    > **Warning:**
    > see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

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

    pvc_name: str = initialize()
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

    input_obj_list, outdir_r_plots = Prepare_Inputs_ATACQC(
        run_flag=NF_Run_Flag, outdir=outdir, aligner=aligner
    )

    input_obj_cov_list, outdir_cov_plots = Prepare_Inputs_Coverages(
        run_flag=NF_Run_Flag, input_dir=outdir, aligner=aligner
    )
    print("Here...")

    map_tasks_cov = map_task(Compress_Coverages_Sample)(IM=input_obj_cov_list)

    map_task_op = map_task(Run_Rscript)(map_input=input_obj_list)
    print(outdir_cov_plots)

    Plotting_Data = Calculate_Plotting_Data(rplots=map_task_op, cov_plots=map_tasks_cov)

    registry_table = UpdateRegistry(d=Plotting_Data, run_name=run_name)

    return outdir
