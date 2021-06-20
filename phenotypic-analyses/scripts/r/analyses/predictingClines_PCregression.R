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

##################################
#### PCA: ENVIRONMENTAL MEANS ####
##################################

# Perform PCA using vegan
pca_enviroMeans <- rda(envMeans_mat, scale = TRUE)
pca_enviroMeans_summary <- summary(pca_enviroMeans)

# Screeplot with broken-stick to examine PC eigenvalues
screeplot(pca_enviroMeans, bstick = TRUE)

# Extract raw scores for first 2 PCs for regression
pca_enviroMeans_rawCityScores <- pca_enviroMeans$CA$u %>% 
  as.data.frame() %>% 
  dplyr::select(PC1, PC2) %>% 
  rename_at(vars(starts_with('PC')), function(x) paste0(x, '_Mean'))


####################################
#### PCA: ENVIRONMENTAL SLOPES #####
####################################

# Perform PCA using vegan
pca_enviroSlopes <- rda(envSlopes_mat, scale = TRUE)
pca_enviroSlopes_summary <- summary(pca_enviroSlopes)

# Screeplot with broken-stick to examine PC eigenvalues
screeplot(pca_enviroSlopes, bstick = TRUE)

# Extract raw scores for first 3 PCs for regression
pca_enviroSlopes_rawCityScores <- pca_enviroSlopes$CA$u %>% 
  as.data.frame() %>% 
  dplyr::select(PC1, PC2, PC3) %>% 
  rename_at(vars(starts_with('PC')), function(x) paste0(x, '_Slopes'))

#######################
#### PC REGRESSION ####
#######################

# Combine PC scores from PCA on Means and Slopes of environmental variables
pca_enviroMeansSlopes_allScores <- bind_cols(pca_enviroMeans_rawCityScores,
                                             pca_enviroSlopes_rawCityScores)

# Create model matrix including interactions of interest
pc_regression_predictors_withInteractions <- model.matrix( ~.^2, data = pca_enviroMeansSlopes_allScores) %>% 
  as.data.frame() %>% 
  dplyr::select(-'(Intercept)') %>% 
  
  # Format interactions
  dplyr::select(-contains(':'),  # Remove all interaction
                ends_with('_Slopes'))

# Create global model
pc_reg_mod <- lm(logOdds_mat ~ ., data = pc_regression_predictors_withInteractions)

# Model selection and averaging using `dredge`
options(na.action = "na.fail")
pc_reg_dredge <- dredge(pc_reg_mod, rank = 'AICc', evaluate = TRUE, extra = c("R^2", "adjR^2"))
options(na.action = "na.omit")

pc_reg_allModels <- as.data.frame(pc_reg_dredge)
pc_reg_topModels <- get.models(pc_reg_dredge, subset = delta < 2)
pc_reg_modelAvg <- model.avg(pc_reg_topModels)
pc_reg_modelAvg_summary <- summary(pc_reg_modelAvg)
pc_reg_topModel_rsquared <- summary(pc_reg_topModels[[1]])$r.squared
