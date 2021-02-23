# Script to generate global table with all data (including collaborator info) for all cities and populations
# 
# Author: James S. Santangelo

# Swansea excluded since too few populations
# Bucaramanga and Santa Marta never sent data
exclude <- c("Swansea", "Bucaramanga", "Santa_Marta")

# Get city centres
allPopMeans <- do.call(rbind, create_df_list(inpath))
city_centres <- allPopMeans %>% 
  dplyr::select(continent, country, city, latitude_city, longitude_city, population) %>% 
  filter(!(city %in% exclude)) %>% 
  distinct()

# Get city stats
city_stats <- read_csv("data/raw/city_data/City_Metrics.csv") %>% 
  dplyr::select(City_Name, Population, Area, Age) %>% 
  rename(city = "City_Name",
         human_population_size = "Population",
         area = "Area",
         age = "Age") %>% 
  mutate(density = human_population_size / area,
         city = as.character(fct_recode(city, "Newhaven" = "New_Haven"))) %>% 
  filter(!(city %in% exclude))

# Cline summary
cline_summary <- read_csv("analysis/supplementary-tables/allCities_HCNslopes_enviroMeansSlopes.csv") %>% 
  dplyr::select(city, betaRLM_freqHCN, pvalRLM_freqHCN, sigRLM) %>% 
  rename("slope_rlm" = "betaRLM_freqHCN",
         "pval_rlm" = "pvalRLM_freqHCN",
         "significant_rlm" = "sigRLM") %>% 
  mutate(slope_rlm = round(slope_rlm, 3)) %>% 
  filter(!(city %in% exclude)) %>% 
  distinct()

# Number of pops and plants. Enviro variables and HCN
enviro_data <- allPopMeans %>% 
  dplyr::select(city, population, total_plants, distance, freqHCN, matches("*Mean$")) %>% 
  dplyr::group_by(city) %>% 
  mutate(num_populations = n(),
         total_num_plants = sum(total_plants),
         transect_length = round(max(distance) - min(distance), 2),
         ID = paste(city, population, sep = "_")) %>% 
  filter(!(city %in% exclude)) %>% 
  ungroup() %>% 
  
  # Below required to remove duplicated rows originating from dataset typos to be fixed
  group_by(city, population) %>% 
  filter(!(n() > 1 & total_plants <= 6)) %>% 
  ungroup() %>% 
  distinct(ID, .keep_all = TRUE)

# Collectors
collectors <- read_csv("data/raw/GLUE_ Team names and affiliations (Responses) - Form Responses 1.csv") %>% 
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

final_table <- city_centres %>% 
  left_join(., city_stats, by = "city") %>% 
  left_join(., cline_summary, by = "city") %>% 
  left_join(., enviro_data, by = c("city", "population")) %>% 
  left_join(., collectors, by = "city") %>% 
  dplyr::select(continent:longitude_city, sampled_by, 
                human_population_size:significant_rlm, freqHCN,
                num_populations:transect_length,
                population, total_plants, contains("Mean")) %>% 
  mutate_at(vars(freqHCN, density, annualAI_Mean, GMIS_Mean:winterNDVI_Mean), round, 3)

write_csv(final_table, path = "analysis/supplementary-tables/globalTable_allInfo_allCities_allPops.csv")
