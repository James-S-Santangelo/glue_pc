# Script identifying which environmental variables predict the strength of clines in HCN
#
# Author: James S. Santangelo and Pedro Peres-Neto

###############
#### SETUP ####
###############

# Need to load a few packages here for cluster execution
library(tidyverse)
library(glmnet)
library(caret)
library(foreach)
source('scripts/r/misc/utilityFunctions.R')

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

#' Estimate Elastic Net model on resampled model matrix
#' 
#' @param i Iteration index
#' @param model_matrix Expanded model matrix with log odds and expanded, standardized predictors
#' 
#' @return Elastic Net model as glmnet object
run_elastic_net_resampled <- function(i, model_matrix){
  set.seed(i)
  elasticNet_model_res <- caret::train(
    y = model_matrix[,1], # Log Odds
    x = model_matrix[,-1], # All predictores
    method = "glmnet",
    metric = "RMSE",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10
  )
  return(elasticNet_model_res)
  
}

#' Extract coefficients for all predictors in elastic net model
#' 
#' @param elasticNet_model Elastic Net model as glmnet object
#' 
#' @return Predictor coefficients as dataframe
extract_coefs_elasticNet <- function(elasticNet_model){
  
  finalModel <- elasticNet_model$finalModel
  lambda <- elasticNet_model$bestTune$lambda
  
  coefs <- as.matrix(coef(finalModel, lambda)) %>% 
    t() %>% 
    as.data.frame()
  
  return(coefs)
}

#' Extract results (e.g., Rsquared, alpha, lambda, etc.) for elastic net model
#' 
#' @param elasticNet_model Elastic Net model as glmnet object
#' 
#' @return Results as data frame
extract_results_elasticNet <- function(elasticNet_model){
  
  # Dataframe with results (RMSE, Rsquared, etc.) for combinations of alpha and lambda
  all_results <- elasticNet_model$results
  
  bestTune_index <- elasticNet_model$bestTune %>% 
    rownames_to_column(., var = 'index') %>% 
    pull(index)
  
  # Get results for best combination of alpha and lambda
  bestTune_results <- all_results[bestTune_index, ]
  
  return(bestTune_results)
}

# Run 'num_reps' elastic net models in parallel for coefficient averaging
registerDoParallel(cores = 24)
num_reps <- 100
elasticNet_list <- foreach(i=1:num_reps, 
                           .verbose = TRUE, 
                           .packages = c('glmnet', 'caret'), 
                           .export = c('model_matrix')) %dopar% 
  run_elastic_net_resampled(i, model_matrix)

# Results and coefficients for observed data
elasticNet_obs_result <- purrr::map_dfr(elasticNet_list, extract_results_elasticNet, .id = 'index')
elasticNet_obs_coefs <- purrr::map_dfr(elasticNet_list, extract_coefs_elasticNet, .id = 'index')

#write_csv(elasticNet_obs_coefs, 'analysis/supplementary-tables/elasticNet_obs_coefs.csv')
#write_csv(elasticNet_obs_result, 'analysis/supplementary-tables/elasticNet_obs_result.csv')
write_csv(elasticNet_obs_coefs, 'analysis/supplementary-tables/elasticNet_obs_coefs2.csv')
write_csv(elasticNet_obs_result, 'analysis/supplementary-tables/elasticNet_obs_result2.csv')

