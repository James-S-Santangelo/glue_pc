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
cline_summary <- read_csv('analysis/supplementary-tables/allCities_logisticReg_coefs.csv') %>% 
  dplyr::select(-continent) %>% 
  mutate(sigLog_Dist = ifelse(pvalLog_Dist < 0.05, "Yes", "No"),
         sigLog_GMIS = ifelse(pvalLog_GMIS < 0.05, "Yes", "No"))

# Number of pops and plants. Enviro variables and HCN
more_city_vars <- allPopMeans %>% 
  dplyr::select(city, population, total_plants, distance, freqHCN) %>% 
  dplyr::group_by(city) %>% 
  mutate(total_plants = ifelse(is.na(total_plants), 20, total_plants)) %>% 
  summarise(num_populations = n(),
            total_num_plants = sum(total_plants),
            transect_length = round(max(distance) - min(distance), 2),
            meanHCN = mean(freqHCN, na.rm = TRUE))

# Collectors
collectors <- read_csv("data/raw/GLUE_collaborators.csv") %>% 
  dplyr::select("City sampled", contains("name"), contains("initial")) %>% 
  rename("city" = "City sampled",
         "Member 1: Middle initial" = "Member 1: Middle initials",
         "Member 2: Middle initial" = "Member 2: Middle initials to be included in publication",
         "Member 3: Middle initial" = "Member 3: Middle initials to be included in publication") %>% 
  rename_at(vars(starts_with("Member")), funs(tolower(.))) %>%  
  rename_at(vars(starts_with("Member")), funs(str_replace_all(., " ", "_"))) %>% 
  rename_at(vars(starts_with("Member")), funs(str_replace_all(., ":", ""))) %>% 
  mutate(city = fct_recode(city, "Portlandme" = "Portland, ME"),
         city = gsub("[;|,|(].*$", "", city),
         city = str_trim(city, side = "right"),
         city = str_replace_all(city, c(" " = "_", "\\." = "",
                                        'ü' = 'u', 'ï' = 'i',
                                        'ë' = 'e', 'ä' = 'a',
                                        'ö' = 'o')),
         city = fct_recode(city, "Frankfurt" = "Frankfurt_am_Main",
                           "Newhaven" = "New_Haven")) %>% 
  filter(!(city %in% exclude),
         !(is.na(city)),
         city != "N/A",
         city != "Sydney") %>% # Sydney excluded since true name of dataset is Parramatta, which is already included
  unite(col = member1, member_1_first_name, member_1_middle_initial, member_1_last_name, sep = " ", remove = TRUE, na.rm = TRUE) %>% 
  unite(col = member2, member_2_first_name, member_2_middle_initial, member_2_last_name, sep = " ", remove = TRUE, na.rm = TRUE) %>% 
  unite(col = member3, member_3_first_name, member_3_middle_initial, member_3_last_name, sep = " ", remove = TRUE, na.rm = TRUE) %>% 
  unite(col = sampled_by, member1, member2, member3, sep = "; ", remove = TRUE, na.rm = TRUE) %>% 
  arrange(city) %>% 
  mutate(sampled_by = gsub("; $", "", sampled_by))

## FINAL TABLE

final_table <- cline_summary %>% 
  left_join(., city_stats, by = "city") %>%
  left_join(., city_centres, by = "city") %>%
  left_join(., more_city_vars, by = "city") %>% 
  left_join(., collectors, by = "city") %>% 
  dplyr::select(continent, Country, city, latitude_city, longitude_city,
                area, pop_size, density, city_age, no_cities, num_populations, 
                total_num_plants, transect_length, betaLog_Dist, betaLog_GMIS, 
                sigLog_Dist, sigLog_GMIS, meanHCN, sampled_by) %>% 
  mutate_if(is.numeric, round, 3)

write_csv(final_table, path = "analysis/supplementary-tables/allCities_stats.csv")
