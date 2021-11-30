## Simple script to get latitudes and longitudes for all sequenced samples

# Load required packages
library(tidyverse)

ss <- read_delim('../genomic-analyses/resources/glue_pc_sampleSheet.txt', delim = '\t') %>%
  rename('population' = 'pop') %>%
  dplyr::select(continent, city, population, individual, site, sample) %>%
  mutate(population = as.character(population))

# Function to load population data
load_pop_data <- function(path){

  df <- read_csv(path) %>%
    dplyr::select(city, country, population, population_latitude, population_longitude) %>%
    mutate(population = as.character(population))
  return(df)
}

inpath <- 'data/clean/popMeans_allCities_withEnviro/'
allPops <- list.files(inpath, full.names = TRUE, pattern = '*.csv') %>%
  map_dfr(., load_pop_data) %>%
  filter(city %in% ss$city)

# Hot Fix for Armidale Sample name clashes
allPops_mod <- allPops %>%
  mutate(population = case_when(city == 'Armidale' ~ str_extract(population, '[0-9]+'),
                                TRUE ~ population))

# Create sample sheet for upload to ENA
ss_forENA <- ss %>%
  left_join(., allPops_mod, by = c('city', 'population')) %>%
  mutate(population_longitude = round(population_longitude, 5),
         population_latitude = round(population_latitude, 5))

write_csv(ss_forENA, '../genomic-analyses/resources/ena_sampleSheet.csv')


