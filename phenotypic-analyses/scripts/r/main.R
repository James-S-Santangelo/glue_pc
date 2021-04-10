# Step-by-step execution of scripts to reproduce results for GLUE's phenotypic data
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
library(interactions)
library(sandwich)
source("scripts/r/utilityFunctions.R")

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
           "analysis/figures/manuscript-panels/",
           "analysis/supplementary-tables/")

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
# (e.g., CGIAR) for each population sampled by collaborators. Scripts are not shown here.
# These scripts were written by Alex Tong. The exception is Impervious surface (GMIS), which 
# we extract using the script below. This requires GMIS raster datasets in a local directory (see script). 
# Datasets can be downloaded from https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1/data-download.
# Accessed April 9, 2021

# Uncomment to execute. Takes about an hour.
# source('scripts/data-extraction/gmis_extraction.R')

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

# Note: Permutations in this script will take a while to run
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

# Estimated increase in HCN frequency per kilometer away from urban core
print(avg_increaseHCN_perKM)

# Anova for random effect of effect of distance varying by city
print(glueClineModel_randEffect_anova)

## Step 4.3: Predicting clines from environment
# What environmental variables drive convergent evolution to cities on a global scale?

set.seed(42)
source("scripts/r/analyses/predictingClines.R")

# Look at residual plots
# Everything looks good
plot(predClines_elasticNet)

# Summary of final Elastiv Net model predicting HCN clines from environmetal data
print(elasticNet_bestTune)  # Alpha and Lambda tuning parameters for Elastic Net. Alpha = 0.9 = Close to full LASSO
print(predClines_elasticNet_summary)

# Anova of final Elastic Net model
print(predClines_elasticNet_anova)

# Summary of Elastic Net model with all Main Effects added back in
print(predClines_elasticNet_withMainEffects_summary)

# Anova of Elastic Net with Main effects back in
print(predClines_elasticNet_withMainEffects_anova)

# Simple slopes analysis for significant winterNDVI_Slope x annualAI_Mean interaction
print(sim_slopes)

# SANITY CHECKS

# Main effect of winterNDVI_Mean goes away when NDSI_Mean is added to model
# likely due to their high correlation (r = 0.93). Let's remove terms with NDSI_Mean
# to see if winterNDVI_Mean comes back. It does.
print(elasticNet_withMainEffects_noNDSI)

# High correlation between winterNDVI_Mean and NDSI_Mean suggests these effects can't be teased apart
# Re-run model selection to see if NDSI_Mean replaces main effect of winterNDVI_Mean. It does.
print(elasticNet_noWinterNDVIMean_summary)

####################################################
#### STEP 5: SUPPLEMENTARY TABLES SNAD DATASETS ####
####################################################

## Step 5.1: Table with best fit cline model summary for each city

source("scripts/r/supplementary-tables/generate_allCities_bestFitModel_clineSummary.R")

## Step 5.2: Table with Means and Robust regression stats for all environmental variables in each city

source("scripts/r/supplementary-tables/generate_enviroMeansSlopes.R")

## Step 5.3: Table combining HCN model stats with environmental means and model stats

source("scripts/r/supplementary-tables/generate_allCities_HCNslopes_enviroMeansSlopes.R")

## Step 5.4: Table with all info (i.e., clines, environmental, collaborator) for all popuations in all cities

source("scripts/r/supplementary-tables/generate_globalTable_allCities_allInfo_allPops.R")
