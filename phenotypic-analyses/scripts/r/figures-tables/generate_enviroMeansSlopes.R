# Script to calculate means and slopes of environmental variables
#
# Author: James S. Santangelo

# Create list with paths to dataframes
inpath <- "data/clean/popMeans_allCities_withEnviro/"
df_list <- create_df_list(inpath)

###################################
#### MEAN ENVIRONMENTAL VALUES ####
###################################

# Calculate means of environmental variables for all cities.
cityEnviroMeans <- purrr::map_df(df_list, calculate_city_eviro_means) %>% 
  rename_at(vars(-city), ~paste0("Mean_", .)) 

#####################################
#### SLOPES ENVIRONMENTAL VALUES ####
#####################################

## Estimate slope of all environmental variables across transect using robust regression

# Vector of environmental variables
envrio_Vars <- c("annualAI_Mean", "annualPET_Mean", "DEM_Mean", "GMIS_Mean", "summerLST_Mean", 
                 "summerNDVI_Mean", "winterLST_Mean", "winterNDVI_Mean", "NDSI_Mean")

# Cross environmental variable vector with population-mean dataframe list for mapping
df_forMapping <- crossing(df_list, envrio_Vars)

# Map robust regression estimation across environmental variables for each city and 'spread'
#   resulting dataframe to get single row per city. 
cityEnviroSlopesRLM <- purrr::map2_dfr(df_forMapping$df_list, df_forMapping$envrio_Vars, rlmStats) %>% 
  pivot_wider(names_from = var, values_from = c(betaRLM, pvalRLM)) %>% 
  rename_at(.vars = vars(ends_with("_Mean")),
            .funs = funs(sub("[_]Mean$", "", .)))

################################
#### MERGE MEANS AND SLOPES ####
################################

enviroMeansSlopes_merged <- cityEnviroMeans %>%
  left_join(., cityEnviroSlopesRLM, by = "city") %>%
  mutate_if(is.numeric, round, 3)

write_csv(enviroMeansSlopes_merged, "analysis/tables/eniroMeansSlopes.csv", col_names = TRUE)

