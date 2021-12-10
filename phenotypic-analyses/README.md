## Analyses of environmental and phenotypic data for GLUE

### Description of repository

This repository contains the analyses of environmental and phenotypic data as part of the GLUE project.

Before running the analyses, the individual-plant phenotype data submitted by collaborators need to be cleaned and standardized. 
This is done using a series of [Python scripts](./scripts/python), which are provided for transparency. However, we have 
distributed the cleaned and standardized datasets with the manuscript so this step is already done. 

### Using the repository

The analyses can be run by following these steps:

1. Open RStudio, navitate to `File` > `Open Project`, and open the `.Rproj` file in this repository.

	- This will load the R project into Rstudio and install the `renv` package if it is not already installed. `renv` is used to manage R pacakge dependencies as part of this project.

2. Run `renv::restore()` to install required packages
 
3. Run [main.R](./scripts/r/main.R). This will create any necessary directories, run analyses, and generate tables and figures. 
This script just calls other scripts that are doing the actual work. Feel free to navigate through these scripts to get a sense of what
they are doing. They should be sufficiently documented to provide an overview of functionality. Briefly, here are the steps in the analysis
pipeline:
    1. Clean latitude and longitude coordinates for city-centers and generate population-mean HCN datasets by city
    2. Clean raw environmental data for each city and merge with population-mean datasets
    3. Analyze environmental data. 
    4. Analyze phenotypic data (i.e., cline analyses)
    5. Predict clines using environmental data
    6. Generate tables and figures. 
