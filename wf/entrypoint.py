from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

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
            "storage_gib": 100,
        }
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: LatchFile, seq_center: typing.Optional[str], read_length: typing.Optional[int], with_control: typing.Optional[bool], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], gtf: typing.Optional[LatchFile], gff: typing.Optional[LatchFile], bwa_index: typing.Optional[str], bowtie2_index: typing.Optional[str], chromap_index: typing.Optional[str], star_index: typing.Optional[str], gene_bed: typing.Optional[LatchFile], tss_bed: typing.Optional[LatchFile], macs_gsize: typing.Optional[float], blacklist: typing.Optional[str], mito_name: typing.Optional[str], save_reference: typing.Optional[bool], ataqv_mito_reference: typing.Optional[str], clip_r1: typing.Optional[int], clip_r2: typing.Optional[int], three_prime_clip_r1: typing.Optional[int], three_prime_clip_r2: typing.Optional[int], trim_nextseq: typing.Optional[int], skip_trimming: typing.Optional[bool], save_trimmed: typing.Optional[bool], keep_dups: typing.Optional[bool], keep_multi_map: typing.Optional[bool], bwa_min_score: typing.Optional[int], save_align_intermeds: typing.Optional[bool], save_unaligned: typing.Optional[bool], narrow_peak: typing.Optional[bool], macs_fdr: typing.Optional[float], macs_pvalue: typing.Optional[float], save_macs_pileup: typing.Optional[bool], skip_peak_qc: typing.Optional[bool], skip_peak_annotation: typing.Optional[bool], skip_consensus_peaks: typing.Optional[bool], skip_deseq2_qc: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_picard_metrics: typing.Optional[bool], skip_plot_profile: typing.Optional[bool], skip_plot_fingerprint: typing.Optional[bool], skip_igv: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], multiqc_methods_description: typing.Optional[LatchFile], fragment_size: typing.Optional[int], keep_mito: typing.Optional[bool], min_trimmed_reads: typing.Optional[int], aligner: typing.Optional[str], skip_merge_replicates: typing.Optional[bool], broad_cutoff: typing.Optional[float], min_reps_consensus: typing.Optional[int], deseq2_vst: typing.Optional[bool], skip_preseq: typing.Optional[bool], skip_ataqv: typing.Optional[bool]) -> None:
    try:
        shared_dir = Path("/nf-workdir")



        ignore_list = [
            "latch",
            ".latch",
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
                *get_flag('input', input),
                *get_flag('fragment_size', fragment_size),
                *get_flag('seq_center', seq_center),
                *get_flag('read_length', read_length),
                *get_flag('with_control', with_control),
                *get_flag('outdir', outdir),
                *get_flag('email', email),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('genome', genome),
                *get_flag('fasta', fasta),
                *get_flag('gtf', gtf),
                *get_flag('gff', gff),
                *get_flag('bwa_index', bwa_index),
                *get_flag('bowtie2_index', bowtie2_index),
                *get_flag('chromap_index', chromap_index),
                *get_flag('star_index', star_index),
                *get_flag('gene_bed', gene_bed),
                *get_flag('tss_bed', tss_bed),
                *get_flag('macs_gsize', macs_gsize),
                *get_flag('blacklist', blacklist),
                *get_flag('mito_name', mito_name),
                *get_flag('save_reference', save_reference),
                *get_flag('keep_mito', keep_mito),
                *get_flag('ataqv_mito_reference', ataqv_mito_reference),
                *get_flag('clip_r1', clip_r1),
                *get_flag('clip_r2', clip_r2),
                *get_flag('three_prime_clip_r1', three_prime_clip_r1),
                *get_flag('three_prime_clip_r2', three_prime_clip_r2),
                *get_flag('trim_nextseq', trim_nextseq),
                *get_flag('min_trimmed_reads', min_trimmed_reads),
                *get_flag('skip_trimming', skip_trimming),
                *get_flag('save_trimmed', save_trimmed),
                *get_flag('aligner', aligner),
                *get_flag('keep_dups', keep_dups),
                *get_flag('keep_multi_map', keep_multi_map),
                *get_flag('bwa_min_score', bwa_min_score),
                *get_flag('skip_merge_replicates', skip_merge_replicates),
                *get_flag('save_align_intermeds', save_align_intermeds),
                *get_flag('save_unaligned', save_unaligned),
                *get_flag('narrow_peak', narrow_peak),
                *get_flag('broad_cutoff', broad_cutoff),
                *get_flag('macs_fdr', macs_fdr),
                *get_flag('macs_pvalue', macs_pvalue),
                *get_flag('min_reps_consensus', min_reps_consensus),
                *get_flag('save_macs_pileup', save_macs_pileup),
                *get_flag('skip_peak_qc', skip_peak_qc),
                *get_flag('skip_peak_annotation', skip_peak_annotation),
                *get_flag('skip_consensus_peaks', skip_consensus_peaks),
                *get_flag('deseq2_vst', deseq2_vst),
                *get_flag('skip_deseq2_qc', skip_deseq2_qc),
                *get_flag('skip_fastqc', skip_fastqc),
                *get_flag('skip_picard_metrics', skip_picard_metrics),
                *get_flag('skip_preseq', skip_preseq),
                *get_flag('skip_plot_profile', skip_plot_profile),
                *get_flag('skip_plot_fingerprint', skip_plot_fingerprint),
                *get_flag('skip_igv', skip_igv),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('skip_qc', skip_qc),
                *get_flag('skip_ataqv', skip_ataqv),
                *get_flag('multiqc_methods_description', multiqc_methods_description)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
        print(flush=True)

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_atacseq", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def nf_nf_core_atacseq(input: LatchFile, seq_center: typing.Optional[str], read_length: typing.Optional[int], with_control: typing.Optional[bool], outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], gtf: typing.Optional[LatchFile], gff: typing.Optional[LatchFile], bwa_index: typing.Optional[str], bowtie2_index: typing.Optional[str], chromap_index: typing.Optional[str], star_index: typing.Optional[str], gene_bed: typing.Optional[LatchFile], tss_bed: typing.Optional[LatchFile], macs_gsize: typing.Optional[float], blacklist: typing.Optional[str], mito_name: typing.Optional[str], save_reference: typing.Optional[bool], ataqv_mito_reference: typing.Optional[str], clip_r1: typing.Optional[int], clip_r2: typing.Optional[int], three_prime_clip_r1: typing.Optional[int], three_prime_clip_r2: typing.Optional[int], trim_nextseq: typing.Optional[int], skip_trimming: typing.Optional[bool], save_trimmed: typing.Optional[bool], keep_dups: typing.Optional[bool], keep_multi_map: typing.Optional[bool], bwa_min_score: typing.Optional[int], save_align_intermeds: typing.Optional[bool], save_unaligned: typing.Optional[bool], narrow_peak: typing.Optional[bool], macs_fdr: typing.Optional[float], macs_pvalue: typing.Optional[float], save_macs_pileup: typing.Optional[bool], skip_peak_qc: typing.Optional[bool], skip_peak_annotation: typing.Optional[bool], skip_consensus_peaks: typing.Optional[bool], skip_deseq2_qc: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_picard_metrics: typing.Optional[bool], skip_plot_profile: typing.Optional[bool], skip_plot_fingerprint: typing.Optional[bool], skip_igv: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], multiqc_methods_description: typing.Optional[LatchFile], fragment_size: typing.Optional[int] = 200, keep_mito: typing.Optional[bool] = False, min_trimmed_reads: typing.Optional[int] = 10000, aligner: typing.Optional[str] = 'bwa', skip_merge_replicates: typing.Optional[bool] = False, broad_cutoff: typing.Optional[float] = 0.1, min_reps_consensus: typing.Optional[int] = 1, deseq2_vst: typing.Optional[bool] = True, skip_preseq: typing.Optional[bool] = True, skip_ataqv: typing.Optional[bool] = False) -> None:
    """
    nf-core/atacseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, fragment_size=fragment_size, seq_center=seq_center, read_length=read_length, with_control=with_control, outdir=outdir, email=email, multiqc_title=multiqc_title, genome=genome, fasta=fasta, gtf=gtf, gff=gff, bwa_index=bwa_index, bowtie2_index=bowtie2_index, chromap_index=chromap_index, star_index=star_index, gene_bed=gene_bed, tss_bed=tss_bed, macs_gsize=macs_gsize, blacklist=blacklist, mito_name=mito_name, save_reference=save_reference, keep_mito=keep_mito, ataqv_mito_reference=ataqv_mito_reference, clip_r1=clip_r1, clip_r2=clip_r2, three_prime_clip_r1=three_prime_clip_r1, three_prime_clip_r2=three_prime_clip_r2, trim_nextseq=trim_nextseq, min_trimmed_reads=min_trimmed_reads, skip_trimming=skip_trimming, save_trimmed=save_trimmed, aligner=aligner, keep_dups=keep_dups, keep_multi_map=keep_multi_map, bwa_min_score=bwa_min_score, skip_merge_replicates=skip_merge_replicates, save_align_intermeds=save_align_intermeds, save_unaligned=save_unaligned, narrow_peak=narrow_peak, broad_cutoff=broad_cutoff, macs_fdr=macs_fdr, macs_pvalue=macs_pvalue, min_reps_consensus=min_reps_consensus, save_macs_pileup=save_macs_pileup, skip_peak_qc=skip_peak_qc, skip_peak_annotation=skip_peak_annotation, skip_consensus_peaks=skip_consensus_peaks, deseq2_vst=deseq2_vst, skip_deseq2_qc=skip_deseq2_qc, skip_fastqc=skip_fastqc, skip_picard_metrics=skip_picard_metrics, skip_preseq=skip_preseq, skip_plot_profile=skip_plot_profile, skip_plot_fingerprint=skip_plot_fingerprint, skip_igv=skip_igv, skip_multiqc=skip_multiqc, skip_qc=skip_qc, skip_ataqv=skip_ataqv, multiqc_methods_description=multiqc_methods_description)

