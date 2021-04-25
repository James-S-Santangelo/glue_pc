# Script to get number of plants per population for Marc's cities

################
#### SETUP #####
################

# This script will read in the individual-plant-level phenotype data from Johnson _et al._ 2018
# and write a single data-frame with the number of plants per population. These datasets were originally
# provided as population-mean datasets where number of plants per pop hadn't been calculated. Because the
# pipeline had been written around these population-mean datasets, it would be too much work to 
# incorporate the individual-plant-level data into the pipeline. As a workaround, I'll calculate the number
# of plants per population from the raw data an include this dataframe in the Git repo, so that these
# data can be merged with the population-mean dataframe that were originally provided. 

# Load individual-level data
johnson_allPlants <- read_csv('data/raw/Johnson et al_20_city_clover 07.12.16_RAW.csv') %>% 
  dplyr::select(City, Site, 'Feigl-Anger') %>% 
  
  # Rename columns
  rename('hcn_result' = 'Feigl-Anger',
         'city' = 'City',
         'population' = 'Site') %>% 
  
  # Extract population as single number
  mutate(population = str_extract(population, pattern = '\\d+')) %>% 
  
  # Remove row with 'NA' as city
  filter(!(is.na(city))) %>% 
    
  # Remap city name
  mutate(city = case_when(city == 'Ge' ~ 'Guelph',
                          city == 'Fe' ~ 'Fergus',
                          city == 'Gt' ~ 'Georgetown',
                          city == 'Ac' ~ 'Acton',
                          city == 'Sf' ~ 'Stratford',
                          city == 'Lo' ~ 'London',
                          city == 'ST' ~ 'Saint_Thomas',
                          city == 'NT' ~ 'New_Tecumseth',
                          city == 'Brt' ~ 'Brantford',
                          city == 'Wo' ~ 'Woodstock',
                          city == 'Wa' ~ 'Waterloo',
                          city == 'Elm' ~ 'Elmira',
                          city == 'elm' ~ 'Elmira',
                          city == 'Or' ~ 'Orangeville',
                          city == 'Ev' ~ 'Everett',
                          city == 'WS' ~ 'Whitchurch-Stouffville',
                          city == 'Brd' ~ 'Bradford',
                          city == 'PH' ~ 'Port_Hope',
                          city == 'An' ~ 'Angus',
                          city == 'Ba' ~ 'Barrie',
                          city == 'Co' ~ 'Cobourg',
                          city == 'co' ~ 'Cobourg')) 

johnson_numPlants <- johnson_allPlants %>% 
  group_by(city, population) %>% 
  summarise(total_plants = sum(!(is.na(hcn_result))),
            numCyanogenic = sum(hcn_result, na.rm = TRUE))  

outpath <- 'data/raw/johnson_2018_plantsPerPop.csv'
write_csv(johnson_numPlants, path = outpath, col_names = TRUE)
