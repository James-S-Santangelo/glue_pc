# Global urban environmental change drives parallel evolution in a cosmopolitan plant
## Authors: XXX

### Abstract

Urbanization dramatically transforms environments in ways that alter the evolution of life. We examined whether urban environmental change drives parallel evolution by sampling 110,019 white clover plants from 6,169 populations in 160 cities spanning diverse climates. Plants were assayed for hydrogen cyanide—a Mendelian antiherbivore defence that also affects tolerance to abiotic stressors. Urban-rural gradients were associated with the evolution of phenotypic clines for hydrogen cyanide in 47% of cities throughout the world. Variation in the strength of clines among cities was explained by environmental changes in drought stress and vegetation that varied among cities. Sequencing 2074 genomes from 26 cities revealed that parallel clines were best explained by adaptive evolution. Our results demonstrate that ongoing urban environmental change is leading to parallel evolution globally.

### Description of repository

This repository contains code and data necessary to reproduce the manuscript's results. The repo can be clones using the following command:

`git clone https://github.com/James-S-Santangelo/glue_pc.git`

Here is a brief description of each subdirectory. Details and documentation can be found in each subdirectory:

- [genomic-analyses](./genomic-analyses): Contains the pipeline used to generate the genomic results derived from low-coverage (\~1X) whole-genome resequencing of 2,074 white clover plants. Uses Conda + Snakemake for reproducibility and pipeline management. 
- [phenotypic-analyses](./phenotypic-analyses): Contains the code used to generate results from the environmental and phenotypic data (i.e., HCN frequencies). Uses Rproject + `renv` for project and dependency management.

