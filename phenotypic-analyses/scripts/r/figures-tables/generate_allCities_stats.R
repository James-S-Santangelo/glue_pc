# Script to generate global table with all data (including collaborator info) for all cities and populations
# 
# Author: James S. Santangelo

# Swansea excluded since too few populations
# Bucaramanga and Santa Marta never sent data
exclude <- c("Swansea", "Bucaramanga", "Santa_Marta")

# Get city centres
inpath <- 'data/clean/popMeans_allCities_withEnviro/'
allPopMeans <- create_df_list(inpath) %>% 
  map(., std_var_zero_one, var = 'GMIS_Mean') %>% 
  do.call(rbind, .)
city_centres <- read_csv("data/clean/latLong_cityCenters_clean.csv")

# Get city stats
city_stats <- read_csv("data/raw/city_data/City_characteristics.csv") %>% 
  dplyr::select(City, area, pop_size, density, city_age, no_cities) %>% 
  rename('city' = 'City') %>% 
  mutate(city = str_replace(city, ' ', '_')) %>% 
  mutate(city = as.character(fct_recode(city, "Newhaven" = "New_Haven"))) %>% 
  filter(!(city %in% exclude))

# Cline summary
cline_summary <- read_csv('analysis/tables/allCities_logisticReg_coefs.csv') %>% 
  dplyr::select(-continent) %>% 
  mutate(sigLog_Dist = ifelse(pvalLog_Dist < 0.05, "Yes", "No"),
         sigLog_GMIS = ifelse(pvalLog_GMIS < 0.05, "Yes", "No"),
         sigLog_hii = ifelse(pvalLog_hii < 0.05, "Yes", "No"))

# Number of pops and plants. Enviro variables and HCN
more_city_vars <- allPopMeans %>% 
  dplyr::select(city, population, total_plants, distance, freqHCN) %>% 
  dplyr::group_by(city) %>% 
  mutate(total_plants = ifelse(is.na(total_plants), 20, total_plants)) %>% 
  summarise(num_populations = n(),
            total_num_plants = sum(total_plants),
            transect_length = round(max(distance) - min(distance), 2),
            meanHCN = mean(freqHCN, na.rm = TRUE))

## FINAL TABLE

final_table <- cline_summary %>% 
  left_join(., city_stats, by = "city") %>%
  left_join(., city_centres, by = "city") %>%
  left_join(., more_city_vars, by = "city") %>% 
  dplyr::select(continent, Country, city, latitude_city, longitude_city,
                area, pop_size, density, city_age, no_cities, num_populations, 
                total_num_plants, transect_length, betaLog_Dist, betaLog_GMIS, 
                betaLog_hii, sigLog_Dist, sigLog_GMIS, sigLog_hii, 
                meanHCN) %>% 
  mutate_if(is.numeric, round, 3)

write_csv(final_table, path = "analysis/tables/allCities_stats.csv")
