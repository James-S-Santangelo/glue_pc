# Simple script for merging extraction and phenotyping data
#
# Author: James S. Santangelo
# Date: Oct. 7 2019


# Load required packages
library(tidyverse)

# Load file with extraction data
extraction_data <- read_csv("resources/illumina-sequencing/library-preps/02_allPlants_Toronto_extractionData_allQC.csv")

# Load file with Toronto plant phenotype data
tor_allPlants <- read_csv("resources/illumina-sequencing/reference/allPlants_Toronto.csv") %>% 
  select(Population, Plant, HCN_Result, Locus.Li, Locus.Ac)

# Load file with Habitat data
habitat <- read_csv("resources/illumina-sequencing/library-preps/01_torontoPops_toExtract.csv") %>% 
  select(Population, Habitat)

allData <- extraction_data %>% 
  left_join(., tor_allPlants, by = c("Population", "Plant")) %>% 
  left_join(., habitat, by = "Population") 

write_csv(allData, path = "resources/illumina-sequencing/library-preps/03_allPlants_Toronto_allData.csv")


