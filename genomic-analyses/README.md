## Analyses of genomic data for GLUE

### TODO

- Write rules to download raw reads and reference genome from Genbank
- Incorporate raw reads for 20 downsampled Toronto samples.
    + These are currently incorporated by calling their BAMs, which were mapped as part of an ongoing project. I find this messy and harder to reproduce since BAMs will not be distributed. 

### Description of repository

This repository contains code necessary to reproduce the genomic analyses in the GLUE manuscript. Raw reads will be archived on Genbank and made available following publication. The pipeline in this repository will generate all of the genomic results directly from the raw reads and the reference genome. The pipeline uses `Conda`, `Snakemake`, and `Singularity` for workflow management and reproducibility. All Snakefiles and directories are well-documented, but here is a brief overview of the pipeline and directories in this repository:

#### Overview of pipeline

1. QC raw reads with [`FastQC`](https://github.com/s-andrews/FastQC)
2. Trim reads with [`fastp`](https://github.com/OpenGene/fastp) and QC trimmed reads
3. Map reads with [`bwa`](https://github.com/lh3/bwa) and sort, index, and mark duplicates using [`SAMtools`](https://github.com/samtools)
4. QC mapped reads using [`bamtools`](https://github.com/pezmaster31/bamtools), [`Qualimap`](http://qualimap.conesalab.org/), [`bamUtil`](https://github.com/statgen/bamUtil), and [`multiQC`](https://github.com/ewels/MultiQC).
5. Estimate diversity and differentiation using [`ANGSD`](https://github.com/ANGSD/angsd) and [`PCAngsd`](https://github.com/Rosemeis/pcangsd).

#### Overview of directories

- [config](./config): Snakemake configuration files for different clusters and programs (e.g., `multiQC`)
- [resources](./resources): Text files used in pipeline (e.g., sample information, chromosomes, etc.)
- [workflow](./workflow): Main Snakemake workflow with rules, environments, scripts, notebooks, and cluster profiles for running the pipeline on Compute Canada SLURM-based clusters.

### Using the pipeline

This pipeline requires `Conda` and `Singularity`:
    - A minimal installation of `Conda` (i.e., Miniconda) can be installed by following the instructions for your platform [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - Installation of `Singularity` requires Admin privileges, but using `Singularity` to run pre-created containers does not. Installation instructions can be found [here](https://sylabs.io/guides/3.5/admin-guide/installation.html). All `Singularity` containers used in this pipeline are avalaible in [this public reposity](https://singularity-hub.org/collections/5044), though they will be automatically pulled and executed by the pipeline. 

Assuming `Conda` is installed, the this repository's `Conda` environment can be replicated by running the following command:

`conda env create -f environment.yaml -n glue`

This will create a `Conda` environment named _glue_ containing a minimal set of dependencies required to run the pipeline (e.g., Python 3.8.6 and Snakemake 5.26.1).

After activating the environment (`conda activate glue`), the pipeline can be executed from the [workflow](./workflow) directory by running a command that looks something like:

`snakemake --use-conda --use-singularity --singularity-args "--bind <path>" --configfile ../config/<configfile> --notemp -j <cores>`

for local execution. Here, `<path>` is the path on the cluster from which files will be read/written (e.g., `/scratch`), `<configfile>` is one of the configfiles in the [config](./config) directory that needs to be modified to match the paths on your system, and `<cores>` is the number of cores available for executing parallel processes. 

For execution on a SLURM cluster, the pipeline can be executed by running:

`snakemake --profile compute-canada --configfile ../config/<configfile>`

Note that the YAML configfiles in the [compute-canada](./workflow/compute-canada/) directory will likely need to be modified to suit your use-case. 
