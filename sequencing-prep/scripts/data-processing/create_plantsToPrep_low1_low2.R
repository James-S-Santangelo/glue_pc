# Script used to select individuals sequenced as part of LOW-1 and LOW-2
#
# Authors: James S. Santangelo and Sophie Koch

#######################
#### LOW1 and LOW2 ####
#######################

## Initially referred to as "Shallow Sequencing", this sampling scheme
## samples 10 urban and 10 rural individuals from each of 49 cities.
## for low coverage (1X) whole genome sequencing (total N = 980)

## Due to time constraints, the sequencing scheme was modified as follows:
## LOW-1  is made up of 500 individuals spread across 25 cities and are used
## as part of GLUE paper 1. LOW-2 represents 150 individuals from each of an 
## additional 24 cities (3 urban, 3 rural per city). These will be part of a later
## paper examining global demography in white clover.

## I first generate the dataset including all 980 plants, which I then split into
## those that have already been sent for sequencing as part of LOW-1, and then sample
## plants that will be sent for sequencing as part of LOW-2

# Load data with extracted Toronto plants
gluePlants <- read_csv("data/clean/extractions/allExtractions.csv",
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
#' 1b. 10 individuals should be selected to span the 5 urban/rural populations
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
plantsToPrep_low1_low2 <- map_dfr(df_list, sampleShallow) %>%
  
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

### LOW-1 ###

## This sheet was created after samples were already sent for sequencing and is generated
## here for completeness

low1_low2_fisrtLanePrepped <- read_csv("resources/20210105_plantsToPrep_Low1_Low2_firstLanePrepped.csv") %>% 
  dplyr::select(continent, city, pop, individual, plantID, 'Batch/lane', "Date prepped", "Bioruptor_label")

low1_prepped <- low1_low2_fisrtLanePrepped %>% 
  rename('lane' = 'Batch/lane',
         "date_prepped" = 'Date prepped') %>% 
  filter(lane == 1 & !(is.na(date_prepped)))

low_1 <- plantsToPrep_low1_low2 %>% 
  filter(plantID %in% low1_prepped$plantID)

# Number of cities
length(unique(low_1$city))

# Number of plants sampled in each population of each city
sampledPlants_byPop_low1 <- low_1 %>% 
  group_by(city, site, pop) %>% 
  tally()

# Number of plants sampled in each site of each city
sampledPlants_bySite_low1 <- low_1 %>% 
  group_by(city, site) %>% 
  tally()

outpath <- 'data/clean/low1/'
dir.create(outpath)
print(sprintf('Plants for LOW1 sequencing saved to %s', outpath))

# Write plants to sample to disk, changing permissions in the process
write_csv(low_1, paste0(outpath, 'plantsToPrep_low1.csv'))

# Write plants by pop and site to disk
write_csv(sampledPlants_byPop_low1, paste0(outpath, 'samplePlants_byPop_low1.csv'))
write_csv(sampledPlants_bySite_low1, paste0(outpath, 'sampledPlants_bySite_low1.csv'))


### LOW-2 ###

## Here is randomly sample 3 urban and rural plants from each city not included in LOW-1

low2_potentialPlants <- plantsToPrep_low1_low2 %>% 
  filter(!(plantID %in% low_1$plantID))

set.seed(2)
low2_toPrep <- low2_potentialPlants %>% 
  # Randomly sample 1 plant per pop
  group_by(city, site, pop) %>% 
  sample_n(1, replace = FALSE) %>% 
  ungroup() %>% 
  # Randomly sample 3 plants per habitat
  group_by(city, site) %>% 
  sample_n(3, replace = FALSE) %>% 
  left_join(., low1_low2_fisrtLanePrepped %>% dplyr::select(plantID, Bioruptor_label), by = 'plantID')

low2_extraPlants <- low2_potentialPlants %>% 
  filter(!(plantID %in% low2_toPrep$plantID))

outpath <- 'data/clean/low2/'
print(sprintf('Plants for LOW2 sequencing saved to %s', outpath))

write_csv(low2_toPrep, paste0(outpath, 'plantsToPrep_low2.csv'))
write_csv(low2_extraPlants, paste0(outpath, 'extraPlants_low2.csv'))
