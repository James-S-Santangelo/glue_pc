# Script to generate table with data for HCN clines plus environmental data for each city
# 
# Author: James S. Santangelo

# Load environmental dataset
enviroMeansSlopes <- read_csv("analysis/supplementary-tables/eniroMeansSlopes.csv")

# Load in datasets  with summaries of HCN clines and environmental data
hcnClinesSummary <- read_csv("analysis/supplementary-tables/allCities_bestFitModel_clineSummary.csv")

# City centres
city_centers <- read_csv("data/clean/latLong_cityCenters_clean.csv") %>%
  dplyr::select(-Country, -continent)

# Get slopes for linear-only models
df_list <- create_df_list("data/clean/popMeans_allCities/")
linSlopesOnly <- purrr::map_dfr(df_list, linearSlopesOnly) %>% 
  ungroup() %>% 
  mutate(city = replace(city, city == "New_Haven", "Newhaven")) %>% 
  mutate_if(is.numeric, round, 3)

# Get slopes and stats from robust regression
betaHCN_RLM <- purrr::map_dfr(df_list, rlmStats, "freqHCN") %>% 
  pivot_wider(names_from = var, values_from = c(betaRLM, pvalRLM)) %>% 
  mutate(city = replace(city, city == "New_Haven", "Newhaven")) %>% 
  mutate(sigRLM = ifelse(pvalRLM_freqHCN < 0.05, "Yes", "No"))

# Add city characteristics
city_data <- read_csv("data/raw/city_data/City_Metrics.csv") %>% 
  dplyr::select(City_Name, Population, Area, City_Num, Age) %>% 
  rename("city" = "City_Name",
         "city_population" = "Population",
         "city_area" = "Area",
         "city_num" = "City_Num",
         "city_age" = "Age")

# Merge datasets
allCities_HCNslopes_enviroMeansSlope <- 
  left_join(hcnClinesSummary, enviroMeansSlopes, by = "city") %>%
  left_join(., city_centers, by = "city") %>% 
  left_join(., linSlopesOnly, by = "city") %>% 
  left_join(., city_data, by = "city") %>% 
  left_join(betaHCN_RLM, by = "city") %>% 
  dplyr::select(city, interceptLinOnly, betaLinOnly, pvalLinOnly, rSquaredLinOnly, direction, sigLinOnly, betaRLM_freqHCN, pvalRLM_freqHCN, sigRLM, modelOrderBestFit, betaBestFit, pvalBestFit, 
                rSquareBestFit, yIntBestFit, predictedBestFit, everything())

write_csv(allCities_HCNslopes_enviroMeansSlope, 
          "analysis/supplementary-tables/allCities_HCNslopes_enviroMeansSlopes.csv",
          col_names = TRUE)
