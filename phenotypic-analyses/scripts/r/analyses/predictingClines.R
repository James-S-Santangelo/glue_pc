# Script identifying which environmental variables predict the strength of clines in HCN
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
                                                                   matches("*Mean$"))
  df_all_popMeans <- rbind(df_all_popMeans, data)
}

df_all_popMeans_excluded <- df_all_popMeans %>% 
  dplyr::filter(!(city %in% c("St_Albert")))

# Run function that returns matrix of slopes of environmental variables for each city
results_statsMatrices <- calculate.stats(all.data = df_all_popMeans_excluded, 
                                         permute = FALSE, 
                                         number.extreme.sites=2)

# Extract slopes of environmental variables.
envSlopes <- results_statsMatrices$slope.matrix %>% 
  dplyr::select(-city) %>% 
  setNames(names(.) %>% stringr::str_replace("Mean", "Slope")) %>% 
  as.matrix()

# Extract mean value of environmental variables for each city
envMeans <- calculate_city_eviro_means(df_all_popMeans_excluded) %>% 
  dplyr::select(-city) %>% 
  drop_na() %>% 
  setNames(paste0(names(.), "_Mean")) %>% 
  as.matrix()

# Get Slopes of HCN clines from Robust regression
HCNslopes <- slope.freqHCN(all.data = df_all_popMeans_excluded, number.extreme.sites=1)
HCNslopes <- HCNslopes$slopes %>% 
  as.data.frame() %>% 
  dplyr::select(rlm) %>%
  rename("HCNslope" = "rlm") %>% 
  as.matrix()

# Distance vector: Distance between urban and rural multivariate environments for each city
distance_vector <- results_statsMatrices$D.UR

# Create dataframe with HCN slopes and environmental data
df_slopes_enviro <- data.frame(HCNslopes, envSlopes, envMeans)

####################################################
#### CORRELATION AMONG ENVIRONMENTAL PREDICTORS ####
####################################################

# Create correlation matrix
envSlopes_corr <- generate_correlation_df(envSlopes)
envMeans_corr <- generate_correlation_df(envMeans)

# Create the plots
pairs(envMeans, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

pairs(envSlopes, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

################################################
#### ENVIRONMENTAL PREDICTORS OF HCN CLINES ####
################################################

# Create Matrix with predictors (including interactions)
# Do not include interactions involving distance vector
# Do not include predictors with interactions among mean environmental variables
predictors <- cbind(distance_vector, envSlopes, envMeans)
predictors_withInteractions <- model.matrix( ~.^2, data = data.frame(scale(predictors)))[,-1]  # remove intercept
predictors_withInteractions_reduced <- predictors_withInteractions %>% 
  as.data.frame() %>% 
  dplyr::select(-contains(':'),  # Remove al interaction
                contains('Slope:')) %>% 
  as.matrix()

# Combine predictors with HCN slopes
model_matrix <- cbind(HCNslopes, predictors_withInteractions_reduced)

# Run Elastic Net model selection
elasticNet_model <- caret::train(
  y = model_matrix[,1], # HCNslopes
  x = model_matrix[,-1], # All predictores
  method = "glmnet",
  metric = "RMSE",
  trControl = trainControl("cv", number = 10, savePredictions = "all"),
  tuneLength = 10
)

# Full list of results
elasticNet_modelResults <- data.frame(elasticNet_model$results)

# Best tuning parameter
elasticNet_bestTune <- elasticNet_model$bestTune
bestAlpha <- elasticNet_bestTune$alpha
bestLambda <- elasticNet_bestTune$lambda

# Final model
EN_finalModel <- elasticNet_model$finalModel

# Extract nonzero coefficients
coef(EN_finalModel, bestLambda)
predictors_nonZero <- which(coef(EN_finalModel, bestLambda)!=0)[2:length(which(coef(EN_finalModel, bestLambda)!=0))]-1
predictors_finalModel <- predictors_withInteractions_reduced[, predictors_nonZero]

# Create data frame for model. Add back in main effects that are not in model
model_df_elasticNet = as.data.frame(cbind(HCNslopes, predictors_finalModel))

# Run final model 
predClines_elasticNet <- lm(HCNslopes ~ ., data = model_df_elasticNet[-1])
predClines_elasticNet_summary <- summary(predClines_elasticNet)
predClines_elasticNet_anova <- Anova(predClines_elasticNet, type = 3)

# Run final model with main effects back in the model.
# Main effects are required to get predicted effects for interactions (Hierarchy principle)
# Model is written out in full to facilitate plotting the interaction
missing_main_effects <- predictors_withInteractions_reduced %>%
  as.data.frame() %>%
  dplyr::select(annualPET_Slope, summerLST_Mean, summerNDVI_Mean,
                DEM_Slope, NDSI_Mean, summerNDVI_Slope, annualAI_Mean,
                winterNDVI_Slope)
model_df_withMainEffects = as.data.frame(cbind(HCNslopes, missing_main_effects, predictors_finalModel))
predClines_elasticNet_withMainEffects <- lm(HCNslopes ~ distance_vector +
                   NDSI_Mean +
                   summerLST_Mean +
                   summerNDVI_Mean +
                   annualAI_Mean +
                   annualPET_Slope +
                   DEM_Slope +
                   summerNDVI_Slope +
                   winterNDVI_Slope +
                   GMIS_Slope +
                   GMIS_Mean + winterNDVI_Mean +
                   annualPET_Slope:summerLST_Mean +
                   annualPET_Slope:summerNDVI_Mean +
                   DEM_Slope:NDSI_Mean +
                   summerNDVI_Slope:annualAI_Mean +
                   summerNDVI_Slope:summerLST_Mean +
                   winterNDVI_Slope:annualAI_Mean +
                   winterNDVI_Slope:summerLST_Mean,
                 data = model_df_withMainEffects)
predClines_elasticNet_withMainEffects_summary <- summary(predClines_elasticNet_withMainEffects)
predClines_elasticNet_withMainEffects_anova <- summary(predClines_elasticNet_withMainEffects)

# Simple slopes analysis for significant winterNDVI_Slope x annualAI_Mean interaction
sim_slopes <- sim_slopes(predClines_elasticNet_withMainEffects, pred = winterNDVI_Slope, modx = annualAI_Mean)

#######################
#### SANITY CHECKS ####
#######################

### Removing NDSI_Mean

## Main effect of winterNDVI_Mean goes away when NDSI_Mean is added to model
## likely due to their high correlation (r = 0.93). Let's remove terms with NDSI_Mean
## to see if winterNDVI_Mean comes back. It does.

model_df_withMainEffects_noNDSI <- model_df_withMainEffects %>%
  dplyr::select(-contains('NDSI_Mean'))

elasticNet_withMainEffects_noNDSI <- lm(HCNslopes ~ ., data = model_df_withMainEffects_noNDSI[-1])
elasticNet_withMainEffects_noNDSI <- summary(elasticNet_withMainEffects_noNDSI)

### Removing winterNDVI_Mean

## High correlation between winterNDVI_Mean and NDSI_Mean suggests these effects can't be teased apart
## Re-run model selection to see if NDSI_Mean replaces main effect of winterNDVI_Mean. It does.

# Remove winterNDVI_Mean from predictor matrix
predictors_withInteractions_reduced_noWinterNDVIMean <- predictors_withInteractions_reduced %>%
  as.data.frame() %>%
  dplyr::select(-contains('winterNDVI_Mean')) %>%
  as.matrix()
model_matrix_noWinterNDVIMean <- cbind(HCNslopes, predictors_withInteractions_reduced_noWinterNDVIMean)

# Run Elastic Net model selection
elasticNet_model_noWinterNDVI_Mean <- caret::train(
  y = model_matrix_noWinterNDVIMean[,1], # HCNslopes
  x = model_matrix_noWinterNDVIMean[,-1], # All predictores
  method = "glmnet",
  metric = "RMSE",
  trControl = trainControl("cv", number = 10, savePredictions = "all"),
  tuneLength = 10
)

# Full list of results
elasticNet_modelResults_noWinterNDVIMean <- data.frame(elasticNet_model_noWinterNDVI_Mean$results)

# Best tuning parameter
elasticNet_model_noWinterNDVI_Mean_bestTune <- elasticNet_model_noWinterNDVI_Mean$bestTune
bestAlpha_noWinterNDVIMean <- elasticNet_model_noWinterNDVI_Mean_bestTune$alpha
bestLambda_noWinterNDVIMean <- elasticNet_model_noWinterNDVI_Mean_bestTune$lambda

# Final model
EN_finalModel_noWinterNDVI <- elasticNet_model_noWinterNDVI_Mean$finalModel

# Extract nonzero coefficients
coef(EN_finalModel_noWinterNDVI, bestLambda_noWinterNDVIMean)
predictors_nonZero_noWinterNDVIMean <- which(coef(EN_finalModel_noWinterNDVI, bestLambda_noWinterNDVIMean)!=0)
predictors_nonZero_noWinterNDVIMean <- predictors_nonZero_noWinterNDVIMean[2:length(which(coef(EN_finalModel_noWinterNDVI, bestLambda_noWinterNDVIMean)!=0))]-1

predictors_finalModel_noWinterNDVIMean <- predictors_withInteractions_reduced_noWinterNDVIMean[, predictors_nonZero_noWinterNDVIMean]

# Create data frame for model. Add back in main effects that are not in model
model_df_noWinterNDVIMean = as.data.frame(cbind(HCNslopes, predictors_finalModel_noWinterNDVIMean))

# Run final model
elasticNet_noWinterNDVIMean <- lm(HCNslopes ~ ., data = model_df_noWinterNDVIMean[-1])
elasticNet_noWinterNDVIMean_summary <- summary(elasticNet_noWinterNDVIMean)
elasticNet_noWinterNDVIMean_anova <- Anova(elasticNet_noWinterNDVIMean, type = 3)
