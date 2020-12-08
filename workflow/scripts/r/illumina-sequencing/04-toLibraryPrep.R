## This script is used for determining which plants will be sequenced
## for the Toronto Pilot project. I additionally select plants for 
# testing and optimizing of the library prep protocol. 

# Load required packages
library(tidyverse)

#### EXPERIMENTAL PLANTS ####

# Load in data with extracted Toronto plants
torPlants <- read_csv("data/illumina-sequencing/library-preps/03_allPlants_Toronto_allData.csv") %>% 
  filter(!is.na(Qubit_conc)) %>% 
  mutate(plant_id = paste(Population, Plant, sep = "_"))

final_vol <- 50 # Volume required for shearing
final_conc <- 10 # Concentration required for shearing
initial_vol <- 40 # Assume I have 40 uL for all sampled (I should have a bit more)

# Need at least this concentration for samples
initial_conc <- (final_vol * final_conc) / initial_vol

# Number of populations by transect and habitat
torPlants %>% 
  group_by(Transect, Distance,Habitat) %>% 
  tally()


# Number of plants in each population that can be extracted
numGoodPlants_byPop <- torPlants %>% 
  group_by(Habitat, Transect, Population, Distance) %>% 
  filter(Qubit_conc > initial_conc) %>% 
  tally() %>% 
  arrange(Distance, .by_group = TRUE) 

# Populations from which to sequence 9 individuals
popsDeepSample <- numGoodPlants_byPop %>% 
  filter(n >= 9) %>% 
  filter(Population %in% c(7, 83, 97) | # Rural populations, manually selected from numGoodPlants_byPop
           Population %in% c(23, 54, 116) | # Suburban populations, manually selected from numGoodPlants_byPop
           Habitat == "Urban") # Taking all 5 urban populations

# Populations from which to sequence 1 individual
popsShallowSample <- numGoodPlants_byPop %>% 
  filter(!(Population %in% popsDeepSample$Population) & # Only sample single plants from populations not deeply sequenced
           Habitat != "Urban" & # Don't need more urban plants
           !(Population %in% c(55, 118))) # Only need 3 extra suburban populations for shallow sampling 

#' Randomly sample plants from a habitat
#' 
#' @param habitat Habitat. One of 'Urban', 'Rural', 'Suburban'.
#' @param num_plants Number of plants to sample.
#' @param pop_df Dataframe with available populations to sample
#' @return Dataframe with randomly sampled plants as rows.
sample_plants <- function(habitat, num_plants, pop_df){
  plants_to_sample <- torPlants %>% 
    filter(Qubit_conc > initial_conc) %>% 
    filter(Population %in% c(pop_df$Population) &
             Habitat %in% habitat) %>% 
    group_by(Population) %>% 
    sample_n(num_plants, replace = FALSE)
  
  return(plants_to_sample)
}

# Radomly sample 9 plants from 3 rural populations
set.seed(43)
ruralPlantsDeepSample <- sample_plants(habitat = "Rural", 
                                       num_plants = 9,
                                       pop_df = popsDeepSample)

# Randomly sample 7 plants from each of 3 suburban populations
suburbanPlantsDeepSample <- sample_plants(habitat = "Suburban", 
                                          num_plants = 7,
                                          pop_df = popsDeepSample)

# Randomly sample 9 plants from each of 5 urban populations
urbanPlantsDeepSample <- sample_plants(habitat = "Urban",
                                       num_plants = 9,
                                       pop_df = popsDeepSample)

# Randomply sample 1 plant from other populations
plantsShallowSample <- sample_plants(habitat = c("Urban", "Suburban", "Rural"),
                                     num_plants = 1,
                                     pop_df = popsShallowSample)

# All plants to prep for Toronto Pilot
allPlant_toPrep <- bind_rows(plantsShallowSample,
                             urbanPlantsDeepSample,
                             suburbanPlantsDeepSample,
                             ruralPlantsDeepSample) %>% 
  mutate(plant_id = paste(Population, Plant, sep = "_")) %>% 
  arrange(Population, Plant) %>% 
  select(City, Population, Plant, plant_id, Plate, Habitat, Qubit_conc, '260_280') %>% 
  mutate(initial_vol = round((final_vol * final_conc) / Qubit_conc, 1),
         TE_vol = round(final_vol - initial_vol, 1)) %>% 
  arrange(Plate, Population, Plant)

write_csv(allPlant_toPrep, path = "data/illumina-sequencing/library-preps/04_allPlants_toPrep.csv", col_names = TRUE)

