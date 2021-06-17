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

# Get percent variance of first two PCs
pca_enviroMeans_eig <- pca_enviroMeans$CA$eig
pca_enviroMeans_percent_var <- pca_enviroMeans_eig * 100 / sum(pca_enviroMeans_eig)
pca_enviroMeans_eig_PC1_varEx <- round(pca_enviroMeans_percent_var[1], 1)  # Percent variance explained by PC1
pca_enviroMeans_eig_PC2_varEx <- round(pca_enviroMeans_percent_var[2], 1)  # Percent variance explained by PC2

# Calculate contribution of each environmental variable to PC2
# https://stackoverflow.com/questions/50177409/how-to-calculate-species-contribution-percentages-for-vegan-rda-cca-objects
pca_enviroMeans_contrib <- round(100*scores(pca_enviroMeans, display = "species", scaling = 0)[,1]^2, 3)

# Extract RDA1 and PC1 species scores and bind contributions
pca_enviroMeans_vars  <- scores(pca_enviroMeans, display = 'species', choices = c(1, 2), scaling = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'var') %>% 
  mutate(contrib = pca_enviroMeans_contrib)

# Plot
pal <- wes_palette("Darjeeling1", 3, type = "continuous")
pca_enviroMeans_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = pca_enviroMeans_vars, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = pca_enviroMeans_vars,
            aes(x = PC1, y = PC2, label = var,
                hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", pca_enviroMeans_eig_PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", pca_enviroMeans_eig_PC2_varEx)) +
  scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
  # scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  # scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  ng1 + theme(legend.position = "top",
              legend.direction="horizontal",
              # legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text = element_text(size=10)) +
  guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
pca_enviroMeans_variableContrib

screeplot(pca_enviroMeans, bstick = TRUE)

# Extract raw scores for PC regression
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

# Get percent variance of first two PCs
pca_enviroSlopes_eig <- pca_enviroSlopes$CA$eig
pca_enviroSlopes_percent_var <- pca_enviroSlopes_eig * 100 / sum(pca_enviroSlopes_eig)
pca_enviroSlopes_eig_PC1_varEx <- round(pca_enviroSlopes_percent_var[1], 1)  # Percent variance explained by PC1
pca_enviroSlopes_eig_PC2_varEx <- round(pca_enviroSlopes_percent_var[2], 1)  # Percent variance explained by PC2

# Calculate contribution of each environmental variable to PC2
# https://stackoverflow.com/questions/50177409/how-to-calculate-species-contribution-percentages-for-vegan-rda-cca-objects
pca_enviroSlopes_contrib <- round(100*scores(pca_enviroSlopes, display = "species", scaling = 0)[,1]^2, 3)

# Extract RDA1 and PC1 species scores and bind contributions
pca_enviroSlopes_vars  <- scores(pca_enviroSlopes, display = 'species', choices = c(1, 2), scaling = 2) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = 'var') %>% 
  mutate(contrib = pca_enviroSlopes_contrib)

# Plot
pal <- wes_palette("Darjeeling1", 3, type = "continuous")
pca_enviroSlopes_variableContrib <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = pca_enviroSlopes_vars, aes(x = 0, xend = PC1, y=0, yend = PC2, color = contrib), 
               size = 2, arrow = arrow(length = unit(0.02, "npc")), alpha = 1) +
  geom_text(data = pca_enviroSlopes_vars,
            aes(x = PC1, y = PC2, label = var,
                hjust = "inward", vjust =  0.5 * (1 - sign(PC1))),
            color = "black", size = 3.5) + 
  xlab(sprintf("PC1 (%.1f%%)", pca_enviroSlopes_eig_PC1_varEx)) + ylab(sprintf("PC2 (%.1f%%)", pca_enviroSlopes_eig_PC2_varEx)) +
  scale_colour_gradientn(colours = rev(pal), breaks = seq(from = 5, to = 25, by = 5)) +
  # scale_x_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  # scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.25)) +
  ng1 + theme(legend.position = "top",
              legend.direction="horizontal",
              # legend.title = element_blank(),
              legend.key.size = unit(0.5, "cm"),
              legend.spacing.x = unit(0.1, "cm"),
              legend.text = element_text(size=10)) +
  guides(color = guide_colourbar(barwidth = 10, barheight = 0.5))
pca_enviroSlopes_variableContrib

screeplot(pca_enviroSlopes, bstick = TRUE)

# Extract raw scores for PC regression
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
summary(pc_reg_modelAvg)
summary(pc_reg_topModels[[1]])
