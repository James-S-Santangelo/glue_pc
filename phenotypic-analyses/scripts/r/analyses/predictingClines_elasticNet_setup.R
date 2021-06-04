# Setip for Elastic Net model execution
#
# Author: James S. Santangelo and Pedro Peres-Neto

###############
#### SETUP ####
###############

# Create dataframe with population-mean HCN and environmental variables for every city
# Needs to be created using `read.csv()` for Pedro's code to work...
inpath <- "data/clean/popMeans_allCities_withEnviro/"
csv.files <- list.files(path = inpath, pattern="*.csv")
df_all_popMeans <- c()
for (i in 1:length(csv.files)){
  data <- read.csv(paste0(inpath, csv.files[i])) %>% dplyr::select(city, 
                                                                   std_distance, 
                                                                   freqHCN,
                                                                   total_plants,
                                                                   matches("*Mean$"))
  df_all_popMeans <- rbind(df_all_popMeans, data)
}

# Run function that returns matrix of slopes of environmental variables for each city
results_statsMatrices <- calculate.stats(all.data = df_all_popMeans, 
                                         permute = FALSE, 
                                         number.extreme.sites=2)

# Extract slopes of environmental variables.
envSlopes <- results_statsMatrices$slope.matrix %>% 
  setNames(names(.) %>% stringr::str_replace("Mean", "Slope")) 

# Extract mean value of environmental variables for each city
envMeans <- calculate_city_eviro_means(df_all_popMeans) %>% 
  drop_na() %>% 
  rename_if(is.numeric, paste0, "_Mean")

# Get change in log-odds of HCN from binomial regression
# These are extracted from the mixed model
logOdds <- coef(glueClineModel_stdDist)$city %>%
  rownames_to_column(var = 'city') %>% 
  dplyr::select(city, std_distance) %>% 
  rename('betaLog' = 'std_distance') %>% 
  filter(city %in% envMeans$city)

# Distance vector: Distance between urban and rural multivariate environments for each city
distance_vector <- results_statsMatrices$D.UR

# Create dataframe with HCN slopes and environmental data
df_slopes_enviro <- left_join(logOdds, envSlopes, by = 'city') %>% 
  left_join(., envMeans, by = 'city')

# Convert to matrices
logOdds_mat <- df_slopes_enviro %>% dplyr::select(betaLog) %>% as.matrix()
envSlopes_mat <- df_slopes_enviro %>% dplyr::select(contains('_Slope')) %>% as.matrix()
envMeans_mat <- df_slopes_enviro %>% dplyr::select(contains('_Mean')) %>% as.matrix()

#####################
#### ELASTIC NET ####
#####################

# Create scaled predictor model matrix, including 2-way interactions
# Include only interactions between slopes of 2 environmental variables,
# and between means and slopes of 2 variables.
predictors <- cbind(distance_vector, envSlopes_mat, envMeans_mat) %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, scale)

predictors_withInteractions <- model.matrix( ~.^2, data = predictors) %>% 
  as.data.frame() %>% 
  dplyr::select(-'(Intercept)') %>% 
  
  # Format interactions
  dplyr::select(-contains(':'),  # Remove all interaction
                contains('Slope:')) %>% 
  as.matrix()

# Combine predictors with log odds into final model matrix
model_matrix <- cbind(logOdds_mat, predictors_withInteractions)