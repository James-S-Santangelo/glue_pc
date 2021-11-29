This directory contains the main components of the Snakemake pipeline. 

- [compute-canada][./compute-canada]: Python scripts and configfiles for execution of pipeline on SLURM-based compute clusters. This profile is a modified version of [this Snakemake SLURM profile](https://github.com/Snakemake-Profiles/slurm).
- [envs](./envs): Conda environment YAML files used throughout pipeline. See specific YAML files for exact versions of software dependencies.
- [rules](./rules): Independent Snakfiles with rules for particular pieces of the pipeline. These are included separately to avoid having an overly long main Snakefile. 
- [scripts](./scripts): Scripts used in the pipeline. 
