# nf-core/atacseq: Crick (CAMP HPC) Configuration

This pipeline has been successfully configured for use on the CAMP HPC cluster at the [The Francis Crick Institute](https://www.crick.ac.uk/).

To use, run the pipeline with `-profile crick`. This will launch the [`crick.config`](../../conf/crick.config) which has been pre-configured with a setup suitable for the CAMP HPC cluster. Using this profile, Nextflow will download a temporary singularity image with all of the required software before execution of the pipeline.

Before running the pipeline you will need to load Nextflow and Singularity using the environment module system on CAMP. You can do this by issuing the commands below:

```
module purge
module load Nextflow/0.32.0
module load Singularity/2.6.0-foss-2016b

nextflow run nf-core/atacseq -profile crick --genome '<genome ID>' --design '<path to your design file>' --email test.user@crick.ac.uk
```

A local copy of the iGenomes resource has been made available on CAMP so you should be able to run the pipeline against any reference available in the [`igenomes.config`](../../conf/igenomes.config) by simply using the `--genome <GENOME_ID>` parameter. Some of the more exotic genomes may not have been downloaded onto CAMP so have a look in the `igenomes_base` path specified in [`crick.config`](../../conf/crick.config), and if your genome of interest isnt present please contact [BABS](mailto:bioinformatics@crick.ac.uk).

Alternatively, if you are running the pipeline regularly for genomes that arent available in the iGenomes resource, we recommend creating a config file with paths to your reference genome indices (see [`reference genomes documentation`](reference_genomes.md) for instructions).

>NB: You will need an account to use the HPC cluster on CAMP in order to run the pipeline. If in doubt contact IT.
