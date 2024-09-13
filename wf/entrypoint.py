import os
import shutil
import subprocess
import sys
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

import requests
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

sys.stdout.reconfigure(line_buffering=True)

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)


def get_flag_defaults(name: str, val: typing.Any, default_val: Optional[typing.Any]):
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: LatchFile
    replicate: int


class Reference(Enum):
    hg19 = "GRCh37 (Homo Sapiens hg19)"
    hg38 = "GRCh38 (Homo Sapiens hg38)"
    mm10 = "GRCm39 (Mus Musculus)"


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


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(str(run_name))

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
            "version": 2,
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
    input: List[SampleSheet],
    run_name: str,
    seq_center: Optional[str],
    read_length: Optional[int],
    with_control: bool,
    outdir: LatchOutputDir,
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
    save_align_intermeds: bool,
    save_unaligned: bool,
    narrow_peak: bool,
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

    if genome_source == "latch_genome_source":
        fasta = os.path.join(
            "s3://latch-public/nf-core/atacseq/",
            latch_genome.name,
            latch_genome.name + ".fa",
        )
        gtf = os.path.join(
            "s3://latch-public/nf-core/atacseq/",
            latch_genome.name,
            latch_genome.name + ".refGene.gtf",
        )
        if aligner.name == "bowtie2":
            bowtie2_index = os.path.join(
                "s3://latch-public/nf-core/atacseq/", latch_genome.name, "bowtie2"
            )
        elif aligner.name == "bwa":
            bwa_index = os.path.join(
                "s3://latch-public/nf-core/atacseq/", latch_genome.name, "bwa"
            )
        elif aligner.name == "star":
            star_index = os.path.join(
                "s3://latch-public/nf-core/atacseq/", latch_genome.name, "star"
            )

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        "docker",
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
        *get_flag_defaults("genome", genome, None),
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
