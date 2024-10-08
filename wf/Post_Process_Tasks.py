import os
import shutil
import subprocess
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from subprocess import run

from dataclasses_json import dataclass_json
from latch import map_task, medium_task, small_task
from latch.account import Account
from latch.functions.messages import message
from latch.ldata.path import LPath
from latch.registry.project import Project
from latch.registry.table import Table
from latch.resources.tasks import custom_task
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from latch_metadata.Compress_Coverage import Compress_Coverages, pa, pq

# from latch_metadata.parameters import Aligner


def latch_listdir(data_directory):
    """
    Helper function to return the elements of a latch directory as a list.
    Input [LatchDir]: data_directory
    Output [List]: files
    """
    files = []
    for x in data_directory.iterdir():
        files.append(x.path)
    return files


class Aligner(Enum):
    bwa = "bwa"
    bowtie2 = "bowtie2"
    star = "star"
    chromap = "chromap"


class Reference(Enum):
    hg19 = "GRCh37 (Homo Sapiens hg19)"
    hg38 = "GRCh38 (Homo Sapiens hg38)"
    mm10 = "GRCm39 (Mus Musculus)"


@dataclass_json
@dataclass
class InputMap_Compress_Coverages:
    data_path: LatchFile
    outPath: LatchDir
    sample: str


@dataclass_json
@dataclass
class InputMap_ATACQC:
    run_flag: str
    sample: str
    genome: str
    bamfile: LatchFile
    baifile: LatchFile
    outdir: LatchDir


@dataclass
class Registry_Obj:
    sample: str
    frag_file: LatchFile
    feat_alignment: LatchFile
    sat_curves: LatchFile
    cov_parquet: LatchDir


@custom_task(cpu=4, memory=48, storage_gib=50)
def Prepare_Inputs_Coverages(
    run_flag: str,
    input_dir: LatchOutputDir,
    aligner: typing.Optional[Aligner] = Aligner.bwa,
) -> (typing.List[InputMap_Compress_Coverages], LatchDir):
    """
    Upon running the nextflow ATAC Seq workflow, the workflow runs a python function to
    make good quality peak visualizations. The function queries the bigwig file of coverages from each
    sample and stores them as parquet file. This function computes an object for running map_task of the
    Compress_Coverages_Sample task. It creates an object of the class InputMap_Compress_Coverages,
    for every sample in the dataset.

    Inputs:
        run_flag [str]: The run_flag returned by the nextflow task
        input_dir [LatchDir]: The LatchDir containing the results from nextflow output.
        aligner [Aligner]: The aligner used by the workflow.

    Outputs:
        input_objects [List[InputMap_Compress_Coverages]]: Objects of type InputMap_Compress_Coverages for every sample
        outdir_cov_plots [LatchDir]: Latch directory containing the coverage results, used for visualization.
        outputs are written to {workflow_output_directory}/{run_name}/{cov_parquet}/
    """
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
            table = pq.read_table(k.local_path)
            sample = Path(k.local_path).name.replace(".parquet", "")
            print(c, sample)
            d[sample] = table["data"]
        pq.write_table(
            pa.table(d),
            os.path.join(out_dir, c + ".parquet"),
        )
    return LatchDir(out_dir, os.path.join(outPath.remote_directory, out_dir))


@custom_task(cpu=4, memory=48, storage_gib=50)
def Compress_Coverages_Sample(IM: InputMap_Compress_Coverages) -> LatchDir:
    """
    Task to compress the coverages further from the bigwig sample as parquet files.
    This is run as map_task, for every sample in the dataset.

    Inputs:
        IM [str]: [List[InputMap_Compress_Coverages]]: Objects of type InputMap_Compress_Coverages for every sample

    Outputs:
        outdir_cov_sample [LatchDir]: Latch directory containing the coverage results per sample.
    """
    bigWig_file = IM.data_path
    sample = IM.sample
    outPath = IM.outPath

    local_dir = sample

    print(bigWig_file, outPath, sample)
    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)

    Compress_Coverages(bigWig_file.local_path, local_dir, sample)

    return LatchDir(local_dir, os.path.join(outPath.remote_directory, local_dir))


@small_task
def Calculate_Plotting_Data(
    rplots: typing.List[LatchDir], cov_plots: typing.List[LatchDir]
) -> (typing.List[Registry_Obj], typing.List[LatchFile]):
    """
    This task prepares an object of type 'Registry_Obj' for every sample. It takes the outputs from the Run_RScript and
    Compress_Covrages_Sample, and these obejcts are then written into a registry table.

    Inputs:
        rplots: List[LatchDir]: LatchDir Paths from running Run_RScript as a map_task
        cov_plots: List[LatchDir]: LatchDir Paths from running Compress_Covrages_Sample as a map_task

    Outputs:
        obj_list: List[Registry_Obj]: List of objects of type Registry_Obj,
        that is used to write files back to the registry tables.
    """
    try:
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
            flag = False
            for k in plots.keys():
                if "featurealignment_coverage.txt" in k:
                    flag = True
            if flag == False:
                cmd = [
                    "echo",
                    "NucleosomeFree\tmononucleosome",
                    ">",
                    "featurealignment_coverage.txt",
                ]
                subprocess.run(" ".join(cmd), shell=True, check=True)
                print(os.listdir(Path()))

                directory = str(Path(plots["Frag_Sizes.txt"].remote_path).parent)
                directory = directory.replace("latch:/", "latch:///")
                print(directory)

                plots["featurealignment_coverage.txt"] = LatchFile(
                    "featurealignment_coverage.txt",
                    f"{directory}/featurealignment_coverage.txt",
                )
                print(plots)
                print("\n")

            d[sample] = plots

        for f in cov_plots:
            sample = Path(f.remote_path).name
            d[sample]["Cov_Plots"] = f

        obj_list = []
        files = []
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
            files.append(d[k]["featurealignment_coverage.txt"])
    except Exception as e:
        message(
            "warning",
            {
                "title": "Could not run Rscript for calculating ATAC seq QC metrics",
                "body": f"Error: {str(e)}",
            },
        )
        return [], []
    return obj_list, files


@custom_task(cpu=4, memory=250, storage_gib=250)
def Run_Rscript(map_input: InputMap_ATACQC) -> LatchDir:
    """
    The workflow calls the RScript that invokes the ATACSeqQC package and persists the data matrices
    to make plots on Latch Plots

    Inputs:
        map_input [List[InputMap_ATACQC]]: Objects of type InputMap_ATACQC for every sample. This is run as a
        map_task and simultaneously processes multiple samples.

    Outputs:
        input_objects [List[InputMap_ATACQC]]: Objects of type InputMap_ATACQC for every sample
        [LatchDir]: Latch directory containing the data matrices used for visualization per sample.
    """

    run_flag = map_input.run_flag
    sample = map_input.sample
    genome = map_input.genome
    bamfile = map_input.bamfile
    baifile = map_input.baifile
    outdir = map_input.outdir

    local_dir = sample + "/"
    shutil.copy(bamfile.local_path, str(Path().resolve()))
    shutil.copy(baifile.local_path, str(Path().resolve()))

    if not os.path.isdir(local_dir):
        os.mkdir(local_dir)

    try:
        cmd_RunRscript = [
            "mamba",
            "run",
            "-n",
            "atacseqqc",
            "Rscript",
            "/root/latch_metadata/ATACSeqQC_Plots.R",
            Path(bamfile.local_path).name,
            local_dir,
            genome,
        ]
        print(" ".join(cmd_RunRscript))
        subprocess.run(" ".join(cmd_RunRscript), shell=True, check=True)
    except Exception as e:
        message(
            "warning",
            {
                "title": "Could not run Rscript for calculating ATAC seq QC metrics",
                "body": f"Error: {str(e)}",
            },
        )

    return LatchDir(local_dir, os.path.join(outdir.remote_directory, sample))


@small_task
def UpdateRegistry(d: typing.List[Registry_Obj], run_name: str) -> str:
    """
    The workflow finally writes the tables necessary for making plots with the plotting layout
    in the registry. The workflow looks for a project with the name, "ATAC_Seq_Results", if it doesn't exist,
    it creates one and inserts a table with the {run_name}.

    Inputs:
        d [List[Registry_Obj]]: Objects of type Registry_Obj ontained from the Calculate_Plotting_Data function.
        run_name, str

    Outputs:
        returns the project id
    """
    try:
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
        print(d)

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
    except Exception as e:
        message(
            "warning",
            {
                "title": "Could not update registry tables",
                "body": f"Error: {str(e)}",
            },
        )
        return "-1"
    return target_table.id


@small_task
def Prepare_Inputs_ATACQC(
    run_flag: str,
    outdir: LatchOutputDir,
    genome: typing.Optional[Reference] = Reference.hg19,
    aligner: typing.Optional[Aligner] = Aligner.bwa,
) -> (typing.List[InputMap_ATACQC], LatchDir):
    """
    Upon running the nextflow ATAC Seq workflow, the workflow uses an R package ATACSeqQC to perform QC analysis.
    This function computes an object for running map_task of the
    Run_Rscript task. It creates an object of the class class InputMap_ATACQC,
    for every sample in the dataset.

    Inputs:
        run_flag [str]: The run_flag returned by the nextflow task
        outdir [LatchDir]: The LatchDir containing the results from nextflow output.
        aligner [Aligner]: The aligner used by the workflow.

    Outputs:
        input_objects [List[InputMap_ATACQC]]: Objects of type InputMap_ATACQC for every sample
        outdir_cov_plots [LatchDir]: Latch directory containing the data matrices used for visualization.
        outputs are written to {workflow_output_directory}/{run_name}/{cov_parquet}/
    """
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
    genome_str = "custom"
    if genome == Reference.hg19:
        genome_str = "hg19"
    if genome == Reference.hg38:
        genome_str = "hg38"

    for f in files:
        print(f)
        if f.endswith(".mLb.mkD.sorted.bam"):
            bamfile = f
            baifile = f + ".bai"
            assert baifile in files, "Missing bai file for " + baifile
            sample = f.split("/")[-1].replace(".mLb.mkD.sorted.bam", "")
            o = InputMap_ATACQC(
                run_flag, sample, genome_str, bamfile, baifile, outdir_r_plots
            )
            object_list.append(o)
    return object_list, outdir_r_plots
