# Scrict to get HCN phenotype and habitat data for sequenced plants

# Load required packages
library(tidyverse)

# Import dataframe with data for 120 plant that were sequenced
sequenced_plants <- read_csv('resources/illumina-sequencing/library-preps/05_allPlants_toPrep_randomized.csv') %>% 
  select(Population, Plant) %>% 
  mutate(Sample = paste0('s_', Population, '_', Plant))

# Import dataframe with all extracted plants containing habitat and HCN phenotype data
phenotype_data <- read_csv('resources/illumina-sequencing/library-preps/03_allPlants_Toronto_allData.csv') %>% 
  select(Population, Plant, HCN_Result, Locus.Li, Locus.Ac, Habitat) %>% 
  mutate(Sample = paste0('s_', Population, '_', Plant)) %>% 
  filter(Sample %in% sequenced_plants$Sample)

# Merge dataframes
dat_out <- left_join(phenotype_data, sequenced_plants, by = c('Sample', 'Population', 'Plant')) %>% 
  select(Sample, Habitat, Population, Plant, everything())

# Write dataframe as Tab-delimited file
write_delim(dat_out, path = 'resources/sequencedPlants_phenotypesHabitat.txt', delim = '\t')
