# Parallel environmental and evolutionary change in response to urbanization on a global scale
## Authors: XXX

### Abstract

Urbanization is a global phenomenon that dramatically alters the environment, but the consistency of this change across cities, and how it affects the evolution of life is poorly understood. We studied how urbanization across 160 cities on all inhabited continents affects environmental change and evolution of the globally distributed plant white clover. Urbanization caused convergent environmental change - cities in different parts of the world tended to be more similar and less variable than nearby nonurban habitats. This urban environmental change drove parallel evolution in an ecologically important trait of white clover in 35% of cities. Urban evolution was driven by urban-nonurban gradients in winter vegetation and snow cover, but the strength and direction of parallel evolution depended on regional aridity. Our results show that urbanization causes environmental homogenization on a global scale, which can be a potent driver of rapid evolution in cities.

### Description of repository

This repository contains code and data necessary to reproduce the manuscript's results. It's actually made up of three seaprate GitHub repositories, which were actively developed side-by-side throughout the project. Two of these repositories are included in this main repo as Git Subtrees. All repositories can be downloaded by running the following command:

`git clone https://github.com/James-S-Santangelo/glue-paper1.git`

Here is a brief description of each subdirectory. Details and documentation can be found in each subdirectory:

- [genomic-analyses](./genomic-analyses): Contains the pipeline used to generate the genomic results derived from low-coverage (~1X) whole-genome resequencing of 520 white clover plants. Uses Conda + Snakemake for reproducibility and pipeline management. 
- [phenotypic-analyses](./phenotypic-analyses): Contains the code used to generate results from the environmental and phenotypic data (i.e., HCN frequencies). Uses Rproject + `renv` for project and dependency management.

Some scripts use files from other repositories, so repos are best executed in the following order:

    1. phenotypic-analyses
    2. sequencing-prep
    3. genomic-analyses
