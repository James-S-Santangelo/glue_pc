# Script to calculate some descriptive statistics for the manuscript
#
# Author: James S. Santangelo

# Load required packages
library(tidyverse)

# Import global table with descriptive stats per city
# Will be provided as supplementary table
globalTable <- read_csv("analysis/supplementary-tables/globalTable_allInfo_allCities_allPops.csv") %>% 
  mutate(total_plants = if_else(is.na(total_plants), 20, total_plants)) # Assume 20 plants for Marc's cities

# Mean number of plants per population with standard errors
globalTable %>% 
  summarise(mean = mean(total_plants),
            total_pops = n(),
            se = mean / sqrt(total_pops),
            min = min(total_plants), 
            max = max(total_plants))

# Mean number of populations per city
globalTable %>% 
  group_by(city) %>% 
  summarise(num_pops = n()) %>% 
  ungroup() %>% 
  summarise(mean_pops = mean(num_pops),
            num_cities = n(),
            se = mean_pops / sqrt(num_cities))

# Total number of plants
num_plants <- globalTable %>% 
  summarise(num_plants = sum(total_plants)) %>% pull()

# Total number of populations
# Each row in global table is a population
num_populations <- globalTable %>% nrow()

# Total number of cities
num_cities <- globalTable %>% distinct(city) %>% nrow()

# Import dataset with slopes and environmental data for each city
allCities_slopes <- read_csv("analysis/supplementary-tables/allCities_HCNslopes_enviroMeansSlopes.csv")

# Percent significant clines
allCities_slopes %>% 
  group_by(sigLinOnly) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)

# Percent clines by direction
allCities_slopes %>% 
  group_by(direction, sigLinOnly) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)
