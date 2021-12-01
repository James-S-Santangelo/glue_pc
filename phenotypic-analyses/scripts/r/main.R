# Step-by-step execution of scripts to reproduce results for GLUE_PC environmental and phenotypic analyses
#
# Author: James S. Santangelo

# This script assumes that all individual-level plant phenotype datasets (i.e., HCN calls) that were
# submitted by GLUE collaborators have been cleaned of errors and standardized for downstream processing.
# The cleaning and standardizing is accomplished using a series of python scripts (../python)

###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)
library(raster)
library(sp)
library(rworldmap)
library(rworldxtra)
library(MASS)
library(heplots)
library(candisc)
library(ggord)
library(InPosition)
library(factoextra)
library(FactoMineR)
library(ggpubr)
library(gridExtra)
library(vegan)
library(wesanderson)
library(Hmisc)
library(glmnet)
library(caret)
library(lme4)
library(lmerTest)
library(sfsmisc)
library(broom)
library(emmeans)
library(ggeffects)
library(interactions)
library(sandwich)
library(patchwork)
library(foreach)
library(doParallel)
library(MuMIn)
source("scripts/r/misc/utilityFunctions.R")

# Create directories
paths <- c("data/clean/environmental_data/annualAI",
           "data/clean/environmental_data/annualPET",
           "data/clean/environmental_data/DEM",
           "data/clean/environmental_data/GMIS",
           "data/clean/environmental_data/NDSI",
           "data/clean/environmental_data/summerLST",
           "data/clean/environmental_data/summerNDVI",
           "data/clean/environmental_data/winterLST",
           "data/clean/environmental_data/winterNDVI",
           "data/clean/popMeans_allCities/",
           "data/clean/popMeans_allCities_withEnviro/",
           "analysis/figures/cline_biplots/",
           "analysis/figures/environmental_biplots/",
           "analysis/figures/main_text/",
           "analysis/figures/main_text/figure2/",
           "analysis/figures/main_text/figure4/",
           "analysis/figures/supplemental/",
           "analysis/tables/")

purrr::walk(paths, dir.create, recursive = T, showWarnings = T)

# Set seed for random processes (e.g., permutations)
set.seed(42)

#################################
#### STEP 1: DATA PROCESSING ####
#################################

## Step 1.1: Clean latitude/longitude coordinates for city centers

source("scripts/r/data-processing/processLatLong_cityCenters.R")

## Step 1.2: Generate population-mean HCN frequency datasets for each city

source("scripts/r/data-processing/generatePopMeans.R")

############################################
#### STEP 2: EXTRACT ENVIRONMENTAL DATA ####
############################################

# This step used a series of custom Python scripts to interface with ArcMap (v. 10.6.1)
# and extract environmental data from Landsat 7/8 images and publicly curated databases
# (e.g., CGIAR) for each population sampled by collaborators. These scripts are messy and likely 
# difficult to run on other machines. Nonetheless, we provide them for transparency in 
# scripts/python/gis_enviro_data_extraction

# There are two exceptions: Impervious surface (GMIS) and Human Influence Index (HII), which are  
# extracted using the scripts below. We provide the raw data from these analyses to prevent users 
# from having to download the rasters unless they really want to. 

# This requires GMIS raster datasets in a local directory (see script). 
# Datasets can be downloaded from https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1/data-download.
# Accessed April 9, 2021

# The HII raster can be obtained from https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic
# Accessed June 25, 2021

# GMIS. Uncomment to execute. Takes about an hour.
# source('scripts/data-extraction/gmis_extraction.R')

# HII. Uncomment to execute. Only takes a courple minutes
# source('scripts/data-extraction/hii_extraction.

######################################
#### STEP 3: MORE DATA PROCESSING ####
######################################

## Step 3.1: Clean/filter environmental data collected for each city

source("scripts/r/data-processing/cleanEnviroData.R")

## Step 3.2: Add population-mean environmental data to population-mean HCN frequency data.

source("scripts/r/data-processing/popMeans_addEnviroData.R")

 ##########################
#### STEP 4: ANALYSES ####
##########################

## Step 4.1: Environmental analyses
# Does urbanization lead to convergent environmental change in cities throughout the world? (Question 1)

# Note: Permutations in this script will take 20 minutes or so
n_perm = 100  # Number of permutations testing for urban/rural difference in mean multivariate environments
source("scripts/r/analyses/enviroAnalyses.R")

# Summary of urban/rural multivariate environment PCA
print(enviroPCA_summary)
print(enviroPCA_summary$cont)  # % variance explained by PC

# Summary of urban/rural multivariate environment permutation test
print(enviroPCA_permute)

# Box-M test for differences in urban/rural multivariate environmental variances
print(enviroVariance_boxM)

# Levene test for differences in urban/rural multivariate environmental variances
print(enviroVariance_leveneLM)
print(enviroVariance_leveneLM_anova)

# Cannonical descriminant analysis for differences in urban/rural multivariate environmental variances
print(enviroVariance_leveneCAN)

## Step 4.2: Cline analyses
# Does convergent urban environmental change cause parallel evolution in an ecologically important trait?
# This script will produce a few warning message, but these are all normal.

source("scripts/r/analyses/clineAnalyses.R")

# Residuals vs. fitted for global GLUE cline model
plot(glueClineModel_stdDist)
hist(residuals(glueClineModel_stdDist))

# Summary of mixed model testing for effect of distance on global scale
print(glueClineModel_stdDist_summary)

# Anova of mixed model testing for effect of distance on global scale
print(glueClineModel_stdDist_anova)

# Marginal and conditional R-squared for global cline model
print(glueClineModel_stdDist_r_squared)

# Estimated percent change in HCN+ across transect. Averaged across cities and continents
print(pred_percent_change)

# LRT for random effect of effect of distance varying by city
print(glueClineModel_stdDist_slopeRE_test)

# Supplementary analysis using GMIS as predictor instead of Distance. Qualitatively similar
print(glueClineModel_stdGmis_anova)
print(glueClineModel_stdGmis_r_squared)
print(glueClineModel_stdGmis_slopeRE_test)

# Supplementary analysis using HII as predictor instead of Distance. Qualitatively similar
print(glueClineModel_stdhii_anova)
print(glueClineModel_stdhii_r_squared)
print(glueClineModel_stdhii_slopeRE_test)

## Step 4.3: Predicting clines from environment
# What environmental variables drive convergent evolution to cities on a global scale?

# Setup for elastic net execution (i.e., defining model matrix)
source("scripts/r/analyses/predictingClines_elasticNet_setup.R")

# Execute Elastic Net. This runs 100 Elastic Net models with 100 different random seeds
# Coefficients were avereaged across the 100 independent runs
# This was executed on a cluster with 24 cores. Feel free to Uncomment and change number of cores for 
# local execution. However, we provide the output of this script as CSV files, which are loaded
# by the elasticNet summary script
num_cores <- 24
# source("scripts/r/analyses/predictingClines_elasticNet_execution.R")

# Summarise output of elastic net model selection and averaging procedure
source("scripts/r/analyses/predictingClines_elasticNet_summary.R")

# Print Elastic Net coefficient summary
print(elasticNet_coefSummary)

# Perform PC regression using environmental data
# Takes a few minutes to run due to `dredge()` model selection and averaging
source('scripts/r/analyses/predictingClines_PCregression.R')

# Summaries of environmental variable PCAs
print(pca_enviroMeans_summary)
print(pca_enviroSlopes_summary)

# Print model averaged coefficeints and R-squared of top model
print(pc_reg_modelAvg_summary)
print(pc_reg_topModel_rsquared)

# Step 4.4: Predict clines from city characteristics
source('scripts/r/analyses/predictingClines_cityCharacteristics.R')

print(finalModel_cityStats)
plot(modlog_popout)
hist(residuals(modlog_popout))

#####################################
#### STEP 5: FIGURES AND TABLES  ####
#####################################

## Step 5.1: Table with best fit cline model summary for each city
source("scripts/r/figures-tables/generate_allCities_logisticReg_coefs.R")

## Step 5.2: Table with Means and Robust regression stats for all environmental variables in each city
source("scripts/r/figures-tables/generate_enviroMeansSlopes.R")

## Step 5.3: Table with city stats, collaborators, HCN slopes
source("scripts/r/figures-tables/generate_allCities_stats.R")

## Step 5.4: Main text figures and tables
## Note: This script uses objects generated in previous scripts
source("scripts/r/figures-tables/main_text_figures.R")

## Step 5.5: Supplemental figures
## Note: This script uses objects generated in previous scripts
source("scripts/r/figures-tables/supplemental_figures.R")

## Step 5.5: Descriptive statistics (e.g., total number of plant, populations, etc.)
source('scripts/r/misc/descriptive_stats.R')

# Mean number of populations per city
print(mean_populations)

# Mean number of plants per population
print(mean_plants_per_pop)

# Total number of plants
print(num_plants)

# Percent significant clines. 
print(percent_sig_clines_logReg)

# Percent significant clines by direction
print(percent_sig_clines_logReg_byDirection)

