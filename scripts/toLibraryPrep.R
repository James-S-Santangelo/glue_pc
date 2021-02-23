# Script to sample individuals for use in GLUE sequencing
#
# Authors: James S. Santangelo and Sophie Koch

###############
#### SETUP ####
###############

# Load required packages
library(tidyverse)

#### SHALLOW SEQUENCING ####

# Shallow sequencing refers to the number of individuals sequenced per
# city, not to sequencing depth In this sampling scheme, 10 urban and 10
# rural individuals are selected and sequenced to 1X coverage from each of
# 49 cities. These individuals will be used for assessing the demographic
# consequences of urbanisation and changes in patterns of genome-wide
# diversity in cities. May also be used for examining selection around the
# HCN loci.

# Load data with extracted Toronto plants
gluePlants <- read_csv("data-clean/allExtractions.csv",
                       col_types = "ccciccnnncccc") %>% 
  # filter(city != "Canberra") %>% 
  mutate(max_qubit = pmax(qubit_1, qubit_2, qubit_3, na.rm = TRUE)) %>% 
  mutate(is_good = ifelse(max_qubit >= 10, 1, 0))

goodPlants <- gluePlants %>% 
  filter(is_good == 1)

# Add number of useable plants by population as column "count"
numToSample_by_pop <- goodPlants %>% 
  group_by(city, site, pop) %>% 
  summarise(numGood = sum(is_good))

#' Select plants to be used for shallow sequencing 
#' 
#' @description 
#' Randomly selects 10 urban and 10 rural individuals from
#' each city for which we have DNA extracted. Individuals are
#' selected as follows:
#' 1a. Randomly select 10 plants with at least 10 ng/uL concentration.
#' 1b. 10 individuals should b selected to span the 5 urban/rural populations
#' (i.e., 2 per population), wherever possible.
#' 2. If fewer than 10 are available, select remaining number of plants
#' with highest concentrations
#' 
#' @param df Dataframe with all extracted plants for a given City and Site
#'   (e.g., Albuquerque, Urban)
#'   
#' @return Dataframe with randomly selected plants for the given City and Site.
sampleShallow <- function(df){
  
  # Get city and site from dataset
  city_name <- unique(df$city)
  site_name <- unique(df$site)
  
  # Assess how many plants meet 10 ng/uL threshold
  df_numPlants <- numToSample_by_pop %>% filter(city == city_name & site == site_name)
  
  # Plants with at least 10 ng/uL concentrations from city and site
  availablePlants_10ng <- goodPlants %>% 
    filter(city == city_name & site == site_name)
  
  total_plants_required <- 10 # Need 10 plants per habitat and city
  total_plants_available <- sum(df_numPlants$numGood) # Number of plants available with 10 ng/uL

  # If there are fewer than 10 total plants for the city/habitat, sample all
  if(total_plants_available < total_plants_required){

    # Sample all available plants with >= 10 ng/uL
    sampledPlants <- availablePlants_10ng %>% 
      sample_n(total_plants_available)
   
    # Calculate number of missing plants
    total_plants_remaining <- total_plants_required - sum(sampledPlants$is_good)

    # Among plants that don't meet 10 ng/uL for current city and site
    # Arrange by decreasing concentration
    # Select top however many plants remaining
    sampledPlants_lessThan10 <- gluePlants %>%
      filter(city == city_name & site == site_name) %>%
      filter(!(plantID %in% sampledPlants$plantID)) %>% # Make sure plant is not already sampled
      arrange(desc(max_qubit)) %>%
      slice(1:total_plants_remaining)
    
    sampledPlants <- bind_rows(sampledPlants, sampledPlants_lessThan10)
    
  # Otherwise group by pop and sample 2 plants from each pop
  }else{
    
    # Sample 2 plants, if possible, otherwise sample 1
    sampledPlants <- availablePlants_10ng %>% 
      group_by(pop) %>% 
      sample_n(ifelse(n() >= 2, 2, 1))
    
    # Calc total number of plants sampled thus far.
    total_plants_sampled <- sum(sampledPlants$is_good)
    
    # If fewer than 10 plants sampled, sample remaining from other pops in city/habitat
    if(total_plants_sampled < total_plants_required){
      
      total_plants_remaining <- total_plants_required - total_plants_sampled
      
      additionalPlants_10ng <- availablePlants_10ng %>% 
        filter(!(plantID %in% sampledPlants$plantID)) %>% 
        sample_n(total_plants_remaining)

      sampledPlants <- bind_rows(sampledPlants, additionalPlants_10ng)
      
    }else{
      sampledPlants <- sampledPlants
    }
      
  }
  return(sampledPlants)
}

# Split all GLUE plant extractions by City and Site, with each
# dataframe as separate list element.
df_list <- goodPlants %>% 
  group_split(city, site)


final_vol <- 25 # Volume required for shearing
final_conc <- 10 # Concentration required for shearing

set.seed(1)

# Select plants from each site and city, combine into single dataframe
plantsToPrep_shallowSample <- map_dfr(df_list, sampleShallow) %>%
  
  # Figure out volume of DNA to remove from stock. 25 uL final volume
  mutate(initial_vol = case_when(
    max_qubit >= 10 ~ round((final_vol * final_conc) / max_qubit, 2),
    
    # If concentration is less than 10, remove 25 uL from stock
    max_qubit < 10 ~ final_vol),
    
    # Volume of TE is 25 uL minus the initial volume removed from stock
         TE_vol = round(final_vol - initial_vol, 2)) %>% 
  
  # Determine whether another library can be prepared in case of failure.
  mutate(leftover = ifelse(initial_vol * 2 < 40, "Yes", "No")) %>% 
  select(-is_good)

# Number of plants sampled in each population of each city
sampledPlants_byPop <- plantsToPrep_shallowSample %>% 
  group_by(city, site, pop) %>% 
  tally()

# Number of plants sampled in each site of each city
sampledPlants_bySite <- plantsToPrep_shallowSample %>% 
  group_by(city, site) %>% 
  tally()

# Write plants to sample to disk, changing permissions in the process
# system("chmod a+w ./data-clean/plantsToPrep_shallowSample.csv")
write_csv(plantsToPrep_shallowSample, "data-clean/plantsToPrep_shallowSample.csv")
# system("chmod a-w ./data-clean/plantsToPrep_shallowSample.csv")

# Write plants by pop and site to disk
write_csv(sampledPlants_byPop, "data-clean/sampledPlants_byPop.csv")
write_csv(sampledPlants_bySite, "data-clean/sampledPlants_bySite.csv")

#### DEEP SEQUENCING ####

# Deep sequencing refers to sampling nd sequencing an additional set of
# individuals from a subset of 26 cities already included in the shallow
# sequencing. For these 26 cities, we will sample and sequence as many
# urban and rural individuals as possible (up to 96 total —— 48 urban and
# 48 rural). These additional individuals will be used for detecting
# genome-wide signatures of selection across the 26 cities.

# Load dataframe with selection of cities to use for Deep sequencing These
# cities were manually selected to balance the presence of clines,
# continental distribution, and number of useable plants. The cities were
# chosen during a meeting between James, Marc, and Rob on February 4, 2020
citySelection <- read_csv("data-raw/deepSample_citySelection.csv") %>% 
  select(city, deep) %>% 
  filter(deep == str_to_lower(deep))

# Select plants for deep sequencing
plantsToPrep_deepSample <- gluePlants %>% 
  
  # Subset to include only cities selected for deep sequencing
  filter(city %in% citySelection$city) %>%
  
  # Remove plants already inluded in shallow sample
  
  filter(!(plantID %in% plantsToPrep_shallowSample$plantID)) %>% 
  
  # Remove plants with no DNA concentration
  filter(max_qubit > 0) %>% 
  
  # Figure out volume of DNA to remove from stock. 25 uL final volume
  mutate(initial_vol = case_when(
    max_qubit >= 10 ~ round((final_vol * final_conc) / max_qubit, 2),
    
    # If concentration is less than 10, remove 25 uL from stock
    max_qubit < 10 ~ final_vol),
    
    # Volume of TE is 25 uL minus the initial volume removed from stock
    TE_vol = round(final_vol - initial_vol, 2)) %>% 
  
  # Determine whether another library can be prepared in case of failure.
  mutate(leftover = ifelse(initial_vol * 2 < 40, "Yes", "No")) %>% 
  select(-is_good)

# Number of plants sampled in each population of each city
sampledPlants_byPop_deep <- plantsToPrep_deepSample %>% 
  group_by(city, site, pop) %>% 
  tally()

# Number of plants sampled in each site of each city
sampledPlants_bySite_deep <- plantsToPrep_deepSample %>% 
  group_by(city, site) %>% 
  tally()

# Write plants to sample to disk, changing permissions in the process
# system("chmod a+w ./data-clean/plantsToPrep_deepSample.csv")
write_csv(plantsToPrep_deepSample, "data-clean/plantsToPrep_deepSample.csv")
# system("chmod a-w ./data-clean/plantsToPrep_deepSample.csv")

# Write plants by pop and site to disk
write_csv(sampledPlants_byPop_deep, "data-clean/sampledPlants_byPop_deepSample.csv")
write_csv(sampledPlants_bySite_deep, "data-clean/sampledPlants_bySite_deepSample.csv")

#### COMPARING SAMPLING VERSIONS ####

compareVersions <- function(df1, df2){
  
  df1 <- df1 %>% 
    ungroup() %>% 
    select(city, plantID)
  
  df2 <- df2 %>% 
    ungroup() %>% 
    select(city, plantID)
  
  out <- anti_join(df1, df2)
  
  return(out)
}

# Compare shallow sample plants selected from scripts to most recent version on Google Drive
mostRecent_shallow <- read_csv("data-clean/previous_shallowSampleSheets/20200530_plantsToPrep_shallowSample.csv")
compareVersions(plantsToPrep_shallowSample, mostRecent_shallow)

# Compare deep sample plants selected from scripts to most recent version on Google Drive
mostRecent_deep <- read_csv("data-clean/previous_deepSampleSheets/20200530_plantsToPrep_deepSample.csv")
compareVersions(plantsToPrep_deepSample, mostRecent_deep)



