# Script to execute Elastic Net model 
#
# Author: James Santangelo

#' Estimate Elastic Net model on resampled model matrix
#' 
#' @param i Iteration index
#' @param model_matrix Expanded model matrix with log odds and expanded, standardized predictors
#' 
#' @return Elastic Net model as glmnet object
run_elastic_net <- function(i, model_matrix){
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
registerDoParallel(cores = num_cores)
num_reps <- 100
elasticNet_list <- foreach(i=1:num_reps, 
                           .verbose = TRUE, 
                           .packages = c('glmnet', 'caret'), 
                           .export = c('model_matrix')) %dopar% 
  run_elastic_net(i, model_matrix)

# Results and coefficients for observed data
elasticNet_obs_result <- purrr::map_dfr(elasticNet_list, extract_results_elasticNet, .id = 'index')
elasticNet_obs_coefs <- purrr::map_dfr(elasticNet_list, extract_coefs_elasticNet, .id = 'index')

write_csv(elasticNet_obs_coefs, 'analysis/tables/elasticNet_obs_coefs.csv')
write_csv(elasticNet_obs_result, 'analysis/ tables/elasticNet_obs_result.csv')


