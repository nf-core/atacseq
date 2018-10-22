# nf-core/atacseq: Crick (CAMP HPC) Configuration

This pipeline has been successfully configured for use on the CAMP HPC cluster at the [The Francis Crick Institute](https://www.crick.ac.uk/).

To use, run the pipeline with `-profile crick`. This will launch the [crick config](../../conf/crick.config) which has been pre-configured with a setup suitable for the CAMP HPC cluster. Using this profile, Nextflow will download a temporary singularity image with all of the required software before execution of the pipeline. **NOTE: You will need an account to use the HPC cluster on CAMP. If in doubt contact IT.**

Before running the pipeline you will need to load Nextflow and Singularity using the environment module system on CAMP. You can do this by issuing the commands below:

```
module use /camp/apps/eb/modules/all
module use /camp/apps/misc/stp/babs/manual/modules/all

module purge
module load nextflow/0.32.0
module load Singularity/2.6.0-foss-2016b
```

A local copy of the iGenomes resource has been made available on CAMP so you should be able to run the pipeline against your reference by simply using the `--genome <GENOME_ID>` parameter. Some of the more exotic genomes may not be available locally so please contact [BABS](bioinformatics@crick.ac.uk) if you would like them to be added.

Alternatively, if you are running the pipeline regularly for custom genomes, we recommend creating a config file with paths to your reference genome indices (see [reference-genomes.md](reference_genomes.md) for instructions).
