# Select plants to use for DEEP-3

################
#### DEEP-3 ####
################

## DEEP-3 refers to a sampling scheme where we sequence more individuals 
## across multiple cities with the goal of searching for parallel signatures
## of selection. These individuals were initially made up of all all individuals
## for which we have DNA extractions from cities that were manually chosen by
## Marc, Rob, and James. 

## Due to time constraints and to maximize the amount of sequencing data, we switched
## the sampling strategy during a meeting in January 2021. We will now sequence all individuals
## for which we have DNA extractions from the same 25 cities included in LOW-1.

# Plants included as part of LOW-1 and LOW-2
low1_prepped <- read_csv('data/clean/low1/plantsToPrep_low1.csv')
low2_prepped <- read_csv('data/clean/low2/plantsToPrep_low2.csv')

# Previous plants chosen as part of DEEP sequencing
deepSample_previousPlants <- read_csv('resources/plantsToPrep_deepSample_withBioruptorTubes.csv')

# Datasheet with DNA extraction data
extraction_data <- read_csv("data/clean/extractions/allExtractions.csv",
                         col_types = "ccciccnnncccc") %>% 
  mutate(max_qubit = pmax(qubit_1, qubit_2, qubit_3, na.rm = TRUE)) %>% 
  mutate(is_good = ifelse(max_qubit >= 10, 1, 0))

goodPlants <- extraction_data %>% 
  filter(is_good == 1)

final_vol <- 25 # Volume required for shearing
final_conc <- 10 # Concentration required for shearing

# Get plants from cities already included in LOW-1
deep3_fromLow1 <- extraction_data %>% 
  # Remove plants from cities not in LOW-1
  filter(city %in% unique(low1_prepped$city)) %>% 
  # Remove plants already sequenced as part of LOW-1
  filter(!(plantID %in% low1_prepped$plantID)) %>% 
  filter(max_qubit > 0)

# Add additional plants from LOW-2 to bring us to ~1900
deep3_fromLow2 <- extraction_data %>% 
  filter(city %in% c("Punta_Arenas", "Palmerston_North", "Warsaw", "Sapporo", "Vancouver")) %>% 
  filter(!(plantID %in% low2_prepped$plantID)) %>% 
  filter(max_qubit > 0)
  

deep3_allPlants <- bind_rows(deep3_fromLow1, deep3_fromLow2) %>% 

  # Figure out volume of DNA to remove from stock. 25 uL final volume
  mutate(initial_vol = case_when(
    max_qubit >= 10 ~ round((final_vol * final_conc) / max_qubit, 2),
  
  # If concentration is less than 10, remove 25 uL from stock
    max_qubit < 10 ~ final_vol),
  
  # Volume of TE is 25 uL minus the initial volume removed from stock
  TE_vol = round(final_vol - initial_vol, 2)) %>% 
  
  # Determine whether another library can be prepared in case of failure.
  mutate(leftover = ifelse(initial_vol * 2 < 40, "Yes", "No")) %>% 
  select(-is_good) %>% 

  arrange(city) %>% 
  
  left_join(., deepSample_previousPlants %>% dplyr::select(plantID, Bioruptor_label), by = 'plantID')

# Get duplicates so they can be corrected
deep3_allPlants %>% 
  filter(duplicated(plantID)) %>% 
  dplyr::select(plantID)

# Add number of useable plants by population as column "count"
numGood_byCitySite <- deep3_allPlants %>% 
  group_by(city, site) %>% 
  summarise(count_deep3Only = n()) %>% 
  mutate(count_withLow1 = count_deep3Only + 10)

outpath <- 'data/clean/deep3/'
print(sprintf('Plants for DEEP3 sequencing save to %s', outpath))
write_csv(deep3_allPlants, paste0(outpath, 'plantsToPrep_deep3.csv'))
write_csv(numGood_byCitySite, paste0(outpath, 'numPlants_byCitySite_deep3.csv'))
