# nf-core/atacseq: Crick (CAMP HPC) Configuration

This pipeline has been successfully configured for use on the CAMP HPC cluster at the [The Francis Crick Institute](https://www.crick.ac.uk/).

To use, run the pipeline with `-profile crick`. This will launch the [`crick.config`](../../conf/crick.config) which has been pre-configured with a setup suitable for the CAMP HPC cluster. Using this profile, Nextflow will download a temporary singularity image with all of the required software before execution of the pipeline.

Before running the pipeline you will need to load Nextflow and Singularity using the environment module system on CAMP. You can do this by issuing the commands below:

```
module purge
module load Nextflow/0.32.0
module load Singularity/2.6.0-foss-2016b

nextflow run nf-core/atacseq -profile crick --genome GRCh37 --design /path/to/design.csv --email test.user@crick.ac.uk
```

A local copy of the iGenomes resource has been made available on CAMP so you should be able to run the pipeline against any reference available in the [`igenomes.config`](../../conf/igenomes.config) by simply using the `--genome <GENOME_ID>` parameter. Some of the more exotic genomes may not have been downloaded onto CAMP so have a look in the `igenomes_base` path specified in [`crick.config`](../../conf/crick.config), and if your genome of interest isnt present please contact [BABS](mailto:bioinformatics@crick.ac.uk).

Alternatively, if you are running the pipeline regularly for genomes that arent available in the iGenomes resource, we recommend creating a config file with paths to your reference genome indices (see [`reference genomes documentation`](reference_genomes.md) for instructions).

If for some reason the pipeline fails to run you can resume it by adding `-resume` to the `nextflow run` command (see [`usage.md`](../usage.md)).  

All of the intermediate files required to run the pipeline will be stored in the `work/` directory. It is recommended to delete this directory after the pipeline has finished successfully because it can get quite large, and all of the main output files will be saved in the `results/` directory anyway. If you wish to keep the `work/` directory for traceability then another solution would be to just delete the largest files. You can execute the commands below from within the `work/` directory:  

```
find ./ -type f -name *.fq.gz -exec rm -rf {} \;
find ./ -type f -name *.sai -exec rm -rf {} \;
find ./ -type f -name *.sam -exec rm -rf {} \;
find ./ -type f -name *.bam -exec rm -rf {} \;
find ./ -type f -name *.bedGraph -exec rm -rf {} \;
find ./ -type f -name *.bigWig -exec rm -rf {} \;
```

>NB: You will not be able to `-resume` the pipeline if you delete the `work/` directory because it contains all of the necessary intermediate files. Please make sure you delete it after the pipeline has finished successfully.

>NB: You will need an account to use the HPC cluster on CAMP in order to run the pipeline. If in doubt contact IT.

>NB: Nextflow will need to submit the jobs via SLURM to the HPC cluster and as such the commands above will have to be executed on one of the login nodes. If in doubt contact IT.
