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
library(doParallel)

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

# Get change in log-odds of HCN from Robust regression
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

####################################################
#### CORRELATION AMONG ENVIRONMENTAL PREDICTORS ####
####################################################

# Create correlation matrix
envSlopes_corr <- generate_correlation_df(envSlopes_mat)
envMeans_corr <- generate_correlation_df(envMeans_mat)

# Create the plots
pairs(envMeans_mat, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

pairs(envSlopes_mat, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

################################################
#### ENVIRONMENTAL PREDICTORS OF HCN CLINES ####
################################################

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

# Define number of reps of elastic net and initialize empty matrix to store coefficients
# num_perms <- 4
# num_reps <- 8
# elasticNet_coefMat <- matrix(0, num_reps, ncol(predictors_withInteractions) + 1)
# colnames(elasticNet_coefMat) <- c('intercept', colnames(predictors_withInteractions))

# Run Elastic Net model with lambda and alpha chosen through repeated (N = 5) 10-fold CV
# Store coefficients for all predictors in matrix after each run. 
set.seed(42)
elasticNet_model <- caret::train(
  y = model_matrix[,1], # Log Odds
  x = model_matrix[,-1], # All predictores
  method = "glmnet",
  metric = "RMSE",
  trControl = trainControl("repeatedcv", repeats = 5, number = 10),
  tuneLength = 10
)

coef(elasticNet_model$finalModel, elasticNet_model$bestTune$lambda) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'predictor') %>% 
  rename('beta' = '1')


run_elastic_net_resampled <- function(i, model_matrix){
  set.seed(i)
  model_matrix_sub <- model_matrix[ sample(nrow(model_matrix), replace = TRUE), ]
  elasticNet_model_res <- caret::train(
    y = model_matrix_sub[,1], # Log Odds
    x = model_matrix_sub[,-1], # All predictores
    method = "glmnet",
    metric = "RMSE",
    trControl = trainControl("repeatedcv", repeats = 5, number = 10),
    tuneLength = 10
  )
  return(elasticNet_model_res)
  
}

num_boot <- 2
registerDoParallel(cores = 2)
ptm <- proc.time()
elasticNet_list <- foreach(i=1:num_boot) %dopar% run_elastic_net_resampled(i, model_matrix)
proc.time() - ptm

extract_coefs_elasticNet <- function(elasticNet_model){
  
  finalModel <- elasticNet_model$finalModel
  lambda <- elasticNet_model$bestTune$lambda
  
  coefs <- as.matrix(coef(finalModel, lambda)) %>% 
    t() %>% 
    as.data.frame()
  
  return(coefs)
  
}

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

elasticNet_obs_result <- extract_results_elasticNet(elasticNet_model) %>% mutate(origin = 'obs')
elasticNet_obs_coefs <- extract_coefs_elasticNet(elasticNet_model) %>% mutate(origin = 'obs')

elasticNet_boot_result <- purrr::map_dfr(elasticNet_list, extract_results_elasticNet, .id = 'index') %>% 
  mutate(origin = 'boot')
elasticNet_boot_coefs <- purrr::map_dfr(elasticNet_list, extract_coefs_elasticNet, .id = 'index') %>% 
  mutate(origin = 'boot')

elasticNet_allResults <- bind_rows(elasticNet_obs_result, elasticNet_boot_result)
elasticNet_allCoefs <- bind_rows(elasticNet_obs_coefs, elasticNet_boot_coefs)

write_csv(elasticNet_allResults, 'analysis/elasticNet_obs_boot_results.csv')
write_csv(elasticNet_allCoefs, 'analysis/elasticNet_obs_boot_coefs.csv')




registerDoParallel(cores = 4)
ptm <- proc.time()
elasticNet_list <- foreach(i=1:num_reps) %dopar% run_elastic_net(i, model_matrix)
proc.time() - ptm



purrr::map_dfr(elasticNet_list, extract_results_elasticNet, .id = 'index')



coefs <- purrr::map_dfr(elasticNet_list, extract_coefs_elasticNet, .id = 'index')

coefs %>% 
  mutate_all(as.numeric) %>% 
  dplyr::select_if(colSums(.) != 0)

tuneLength <- 20
for(i in 1:10) seeds[[i]]<- sample.int(n=1000, tuneLength)
#for the last model
for(i in 1:num_reps){
  set.seed(i)
  print(i)
  elasticNet_model <- caret::train(
    y = model_matrix[,1], # Log Odds
    x = model_matrix[,-1], # All predictores
    method = "glmnet",
    metric = "RMSE",
    trControl = trainControl("cv", number = 10),
    tuneLength = tuneLength
  )
  coefs <- coef(elasticNet_model$finalModel, elasticNet_model$bestTune$lambda)
  elasticNet_coefMat[i,] <- as.matrix(coefs)
}
all.equal(predict(elasticNet_model1, type = 'raw'), predict(elasticNet_model2, type = 'raw'))

coef1 <- coef(elasticNet_model1$finalModel)
coef2 <- coef(elasticNet_model2$finalModel)
elasticNet_model1$bestTune





# Extract only coefficients that are non-zero in at least one model
elasticNet_coefMat_nonZero <- elasticNet_coefMat %>% 
  as.tibble() %>% 
  dplyr::select_if(colSums(.) != 0)

# Take mean of coefficients across 100 models and count number of models in which 
# it was non zero.
elasticNet_coefMeans <- elasticNet_coefMat_nonZero %>% 
  as.tibble() %>% 
  pivot_longer(names(.), names_to = 'predictor', values_to = 'beta') %>% 
  mutate(is_zero = ifelse(beta == 0, 0, 1)) %>% 
  group_by(predictor) %>% 
  summarise(mean_beta = mean(beta),
            count_nonZero = sum(is_zero)) %>% 
  arrange(desc(abs(mean_beta)))

# Dataframe with only final predictors
coeffs_noIntercept_names <- elasticNet_coefMat_nonZero %>% 
  dplyr::select(-intercept) %>% 
  colnames()

finalModel_df <- predictors_withInteractions %>% 
  as.tibble() %>%
  dplyr::select(one_of(coeffs_noIntercept_names)) %>% 
  as.matrix()

# set.seed(42)
# elasticNet_model <- caret::train(
#   y = model_matrix[,1], # Log Odds
#   x = model_matrix[,-1], # All predictores
#   method = "glmnet",
#   metric = "RMSE",
#   trControl = trainControl("cv", number = 10, savePredictions = "all"),
#   tuneLength = 25
# )
# 
# # Full list of results
# elasticNet_modelResults <- data.frame(elasticNet_model$results)
# 
# # Best tuning parameter
# elasticNet_bestTune <- elasticNet_model$bestTune
# bestAlpha <- elasticNet_bestTune$alpha
# bestLambda <- elasticNet_bestTune$lambda
# 
# # Final model
# EN_finalModel <- elasticNet_model$finalModel
# 
# # Extract nonzero coefficients
# EN_final_modelCoefs <- coef(EN_finalModel, bestLambda)
# predictors_nonZero <- which(EN_final_modelCoefs!=0)[2:length(which(EN_final_modelCoefs!=0))]-1
# predictors_finalModel <- predictors_withInteractions[, predictors_nonZero]
# 
# # Create data frame for model. Add back in main effects that are not in model
# model_df_elasticNet = as.data.frame(cbind(logOdds_mat, predictors_finalModel))
# 
# # Run final model 
# predClines_elasticNet <- lm(betaLog ~ ., data = model_df_elasticNet)
# 
# # Diagnostics. A couple outliers
# plot(predClines_elasticNet)
# hist(residuals(predClines_elasticNet))
# 
# # Model summary
# predClines_elasticNet_summary <- summary(predClines_elasticNet)
# predClines_elasticNet_anova <- Anova(predClines_elasticNet, type = 3)

