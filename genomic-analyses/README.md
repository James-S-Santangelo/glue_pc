## Analyses of genomic data for GLUE_PC

### Description of repository

This repository contains code necessary to reproduce the genomic analyses in the GLUE_PC manuscript. BAM files have been deposited on the European Nucleotide Archive (ENA, Study accession #PRJEB48967). The *from_bams* branch in this repository can be used to reproduce the manuscript's results from these BAM files. The _master_ branch runs instead from the raw reads and includes read trimming, mapping, and QC which are largely performed by a [separate pipeline](https://github.com/James-S-Santangelo/glue_dnaSeqQC) that is automatically incorporated here as a Snakemake module. The raw reads are still being used for multiple ongoing projects associated with GLUE but are available upon request. The pipeline uses `Conda`, `Snakemake`, and `Singularity` for workflow management and reproducibility. All Snakefiles and directories are well-documented, but here is a brief overview of the pipeline and directories in this repository:

#### Overview of pipeline

The following steps are perfomed using rules from the separate [dnaSeqQC pipeline](https://github.com/James-S-Santangelo/glue_dnaSeqQC):

1. QC raw reads with [`FastQC`](https://github.com/s-andrews/FastQC)
2. Trim reads with [`fastp`](https://github.com/OpenGene/fastp) and QC trimmed reads
3. Map reads with [`bwa`](https://github.com/lh3/bwa) and sort, index, and mark duplicates using [`SAMtools`](https://github.com/samtools)
4. QC mapped reads using [`bamtools`](https://github.com/pezmaster31/bamtools), [`Qualimap`](http://qualimap.conesalab.org/), [`bamUtil`](https://github.com/statgen/bamUtil), and [`multiQC`](https://github.com/ewels/MultiQC).

The following steps are performed using rules in this repo's Snakemake pipeline. If using the pipeline on the *from_bams* branch, these are the steps that will run.

1. Estimate urban and rural diversity, differentiation, and population structure using [`ANGSD`](https://github.com/ANGSD/angsd) and [`PCAngsd`](https://github.com/Rosemeis/pcangsd).
2. Estimate differentiation of HCN and underlying loci relative to neutral expectations

#### Overview of directories

- [analyses](./analyses): Summary datasets and figures generated as part of this pipeline.
- [config](./config): Snakemake configuration files for different clusters.
- [notebooks](./notebooks): Jupyter Notebooks detailing analyses of diversity, population structure, and HCN/Ac/Li differentiation.
- [resources](./resources): Text files used in pipeline (e.g., sample information, chromosomes, etc.)
- [workflow](./workflow): Main Snakemake workflow with rules, environments, scripts, notebooks, and cluster profiles for running the pipeline on Compute Canada SLURM-based clusters.

### Using the pipeline

This pipeline requires `Conda` and `Singularity`:

- A minimal installation of `Conda` (i.e., Miniconda) can be installed by following the instructions for your platform [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Installation of `Singularity` requires Admin privileges, but using `Singularity` to run pre-created containers does not. Installation instructions can be found [here](https://sylabs.io/guides/3.5/admin-guide/installation.html). All `Singularity` containers used in this pipeline are avalaible in [this public reposity](https://cloud.sylabs.io/library/james-s-santangelo), though they will be automatically pulled and executed by the pipeline. 

Assuming `Conda` is installed, the this repository's `Conda` environment can be replicated by running the following command:

`conda env create -f environment.yaml -n glue_pc`

This will create a `Conda` environment named _glue\_pc_ containing a minimal set of dependencies required to run the pipeline (e.g., Python 3.8.6 and Snakemake 6.9.1). This environment additinally contains dependencies to run Jupter Notebooks (e.g., Jupyter + R and associated packages).

After activating the environment (`conda activate glue_pc`), the pipeline can be executed from the [workflow](./workflow) directory by running a command that looks something like:

`snakemake --use-conda --use-singularity --singularity-args "--bind <path>" --configfile ../config/<configfile> --notemp -j <cores>`

for local execution. Here, `<path>` is the path on the cluster from which files will be read/written (e.g., `/scratch`), `<configfile>` is one of the configfiles in the [config](./config) directory that needs to be modified to match the paths on your system, and `<cores>` is the number of cores available for executing parallel processes. 

For execution on a SLURM cluster, the pipeline can be executed by running:

`snakemake --profile compute-canada --configfile ../config/<configfile>`

Note that the YAML configfiles in the [compute-canada](./workflow/compute-canada/) directory will likely need to be modified to accomodate the paths on your cluster.

If the entire repository has been cloned by running `git clone https://github.com/James-S-Santangelo/glue_pc.git`, the *from_bams* branch can be retrieved by running:

1. `git fetch origin from_bams`
2. `git checkout from_bams`

Because it can take a long time to run the pipeline from the BAM files through to all of the analyses, we've included [summary datasets](./analyses/tables/) for all major analyses. These summary datasets are used in the Jupyter notebooks, and the code that was used to generate them is also provided therein for transparency. 
