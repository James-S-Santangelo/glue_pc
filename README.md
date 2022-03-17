[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5765252.svg)](https://doi.org/10.5281/zenodo.5765252)

## Global Urban Evolution Project – Parallel Clines (GLUE_PC)
### Manuscript: Global urban environmental change drives adaptation in white clover

This manuscript has been published in [Science](https://www.science.org/doi/10.1126/science.abk0989)

### Abstract

Urbanization transforms environments in ways that alter biological evolution. We examined whether urban environmental change drives parallel evolution by sampling 110,019 white clover plants from 6169 populations in 160 cities globally. Plants were assayed for a Mendelian antiherbivore defense that also affects tolerance to abiotic stressors. Urban-rural gradients were associated with the evolution of clines in defense in 47% of cities throughout the world. Variation in the strength of clines was explained by environmental changes in drought stress and vegetation cover that varied among cities. Sequencing 2074 genomes from 26 cities revealed that the evolution of urban-rural clines was best explained by adaptive evolution, but the degree of parallel adaptation varied among cities. Our results demonstrate that urbanization leads to adaptation at a global scale.

### Description of repository

This repository contains code and data necessary to reproduce the manuscript's results. The repo can be clones using the following command:

`git clone https://github.com/James-S-Santangelo/glue_pc.git`

Here is a brief description of each subdirectory. Details and documentation can be found in each subdirectory:

- [genomic-analyses](./genomic-analyses): Contains the pipeline used to generate the genomic results derived from low-coverage (\~1X) whole-genome resequencing of 2,074 white clover plants. Uses Conda + Snakemake for reproducibility and pipeline management. 
- [phenotypic-analyses](./phenotypic-analyses): Contains the code used to generate results from the environmental and phenotypic data (i.e., HCN frequencies). Uses Rproject + `renv` for project and dependency management.

Note: The phenotypic analyses directory produces files used by the Jupyter Notebooks in the genomics analyses directory so it should be run first. However, all required files are present in the GitHub repo and archived data repositories. 

