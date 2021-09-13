# Main script for running all steps of pipeline
#
# Author: James S. Santangelo

###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)

paths <- c('data/clean/extractions/',
           'data/clean/low2/',
           'data/clean/deep3/')

purrr::walk(paths, dir.create, recursive = TRUE, showWarnings = TRUE)

#######################################
#### STEP 1: Clean extraction data ####
#######################################

# Run script to clean extraction data for analysis
source('scripts/data-processing/clean_extractionData.R')

################################
#### STEP 2a: LOW1 and LOW2 ####
################################

# Create datasheets with plants to sequence as part of LOW1 and LOW2
source('scripts/data-processing/create_plantsToPrep_low1_low2.R')

##############################################
#### STEP 2b: LOW1 library concentrations ####
##############################################

# Add Illumina indices and clean post-PCR library concentrations
source('scripts/data-processing/create_low1_libraryConcentrations.R')

####################################
#### STEP 2c: LOW1 sample sheet ####
####################################

# Create sample sheet with habitat info for downstream genomics analyses
source('scripts/data-processing/create_low1_sampleSheet.R')

#######################
#### STEP 3: DEEP3 ####
#######################

# Create datasheets with plants to sequence as part of DEEP3
source('scripts/data-processing/create_plantsToPrep_deep3.R')

#############################################
#### STEP 4: Add lanes to DEEP3 and LOW2 ####
#############################################

# Run script to concatenate and add sequencing lanes to DEEP3 and LOW2
source('scripts/data-processing/concat_deep3_low2_addLanes.R')

#####################################################
#### STEP 5a: DEEP3 LANE1 library concentrations ####
#####################################################

# Add Illumina indices and clean post-PCR library concentrations
# for the first lane of DEEP3 sequencing
source('scripts/data-processing/create_deep3_lane1_libraryConcentrations.R')

############################################
#### STEP 5b: DEEP3, LANE1 sample sheet ####
############################################

# Create sample sheet with habitat info for downstream genomics analyses
source('scripts/data-processing/create_deep3_lane1_sampleSheet.R')

#####################################################
#### STEP 6a: DEEP3 LANE2 library concentrations ####
#####################################################

# Add Illumina indices and clean post-PCR library concentrations
# for the first lane of DEEP3 sequencing
source('scripts/data-processing/create_deep3_lane2_libraryConcentrations.R')

############################################
#### STEP 6b: DEEP3, LANE2 sample sheet ####
############################################

# Create sample sheet with habitat info for downstream genomics analyses
source('scripts/data-processing/create_deep3_lane2_sampleSheet.R')

##########################################################
#### STEP 7a: DEEP3/LOW2 LANE3 library concentrations ####
##########################################################

# Add Illumina indices and clean post-PCR library concentrations
# for the first lane of DEEP3 sequencing
source('scripts/data-processing/create_deep3_low2_lane3_libraryConcentrations.R')

############################################
#### STEP 7b: DEEP3, LANE2 sample sheet ####
############################################

# Create sample sheet with habitat info for downstream genomics analyses
source('scripts/data-processing/create_deep3_low2_lane3_sampleSheet.R')
