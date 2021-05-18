
num_boot <- 1000
elasticNet_obs_boot_coefs <- read_csv('analysis/elasticNet_obs_boot_coefs.csv')
elasticNet_obs_boot_results <- read_csv('analysis/elasticNet_obs_boot_results.csv')

elasticNet_obs_coefs <- elasticNet_obs_boot_coefs %>% 
  filter(origin == 'obs') %>% 
  dplyr::select(-index) %>% 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0)
elasticNet_obs_coefs_names <- colnames(elasticNet_obs_coefs)

elasticNet_boot_coefs <- elasticNet_obs_boot_coefs %>% 
  filter(origin == 'boot') %>% 
  dplyr::select(-index) %>% 
  dplyr::select(one_of(elasticNet_obs_coefs_names))

elasticNet_boot_coefs %>% 
  summarise_each(funs(sum(.==0))) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'predictor') %>% 
  rename('num_zero' = 'V1') %>% 
  mutate(p = num_zero / num_boot)


pairs(elasticNet_boot_coefs %>% dplyr::select(-origin, -`(Intercept)`) %>% as.matrix(), 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
hist(elasticNet_boot_coefs$winterNDVI_Mean)
plot(elasticNet_model$results)

 
library(selectiveInference)
lambda <- 0.245
alpha <- 0.6
results <- glmnet(x = model_matrix[,-1],
                  y = model_matrix[,1], 
                  family = 'gaussian', 
                  alpha = alpha, 
                  lambda = lambda, 
                  intercept = TRUE)
sigma <- estimateSigma(x = model_matrix[,-1],
                       y = model_matrix[,1],
                       intercept = TRUE,
                       standardize = FALSE)
beta <- coef(elasticNet_model$finalModel, x = model_matrix[,-1], s = lambda, y = model_matrix[,1])[-1]
out = fixedLassoInf(x = model_matrix[,-1],
                    y = model_matrix[,1],
                    beta = beta,
                    lambda = lambda,
                    sigma=sigma$sigmahat,
                    alpha = 0.05,
                    tol.kkt = 20)

predict(results, newx = new_model_matrix,
  type = 'response', s = results$lambda)

test <- model_matrix %>% 
  as_tibble() %>% 
  dplyr::select(-betaLog) %>% 
  mutate_at(vars(-(`summerNDVI_Slope:annualAI_Mean`)), mean) %>% 
  mutate("summerNDVI_Slope:annualAI_Mean" = seq(from = min(`summerNDVI_Slope:annualAI_Mean`), 
                    to = max(`summerNDVI_Slope:annualAI_Mean`),
                    length.out = n())) %>% 
  as.matrix()
new_model_matrix <- Matrix(test)

cbind(elasticNet_obs_coefs_names[2:7], out$pv) %>% as_tibble() 

coef(elasticNet_model$finalModel, x = model_matrix[,-1], s = lambda, y = model_matrix[,1])[-1]
coef(elasticNet_model$finalModel, s = lambda)

df <- data.frame(model_matrix)


mod <- lm(betaLog ~ annualPET_Slope + NDSI_Mean + winterNDVI_Mean +
            GMIS_Mean + summerNDVI_Mean + summerNDVI_Slope + annualAI_Mean +
            annualPET_Slope:GMIS_Mean + annualPET_Slope:summerNDVI_Mean + 
            summerNDVI_Slope:annualAI_Mean, data = df)
summary(mod)

# Define low, mean, and high AI categories
LSThigh <- round(mean(df$GMIS_Mean) + sd(df$GMIS_Mean), 1)
LSTmean <- round(mean(df$GMIS_Mean), 1)
LSTlow <- round(mean(df$GMIS_Mean) - sd(df$GMIS_Mean), 1)

# Get dataframe with predicted lines from Elastic Net model for effects of winterNDVI_Slope on HCN slopes for each of 
# 3 levels summerLST_Mean
range(df$summerNDVI_Slope) # Used to parameterize x-axis range in list below. 
LST_wNDVI_list <- list(annualPET_Slope = seq(from = -2, to = 3, by = 0.1), GMIS_Mean = c(LSTlow, LSTmean, LSThigh))
LST_wNDVI_df <- emmip(mod, GMIS_Mean~annualPET_Slope, at = LST_wNDVI_list, CIs=TRUE, plotit=FALSE)

# Color palette
LSTlow_col = wes_palette("Zissou1", 5, type = "discrete")[5]
LSTmean_col = wes_palette("Zissou1", 5, type = "discrete")[3]
LSThigh_col = wes_palette("Zissou1", 5, type = "discrete")[1]
cols = c(LSTlow_col, LSThigh_col)

LST_wNDVI_df$fGMIS_Mean <- factor(LST_wNDVI_df$GMIS_Mean)
LST_wNDVI_plot <- LST_wNDVI_df %>% 
  filter(fGMIS_Mean != '0') %>% 
  ggplot(., aes(x = annualPET_Slope, y = yvar)) + 
  geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = fGMIS_Mean), alpha = 0.25) +
  geom_line(aes(color = fGMIS_Mean), size = 2) +
  ylab("Predicted HCN slope") + xlab("Slope of annual PET") +
  scale_color_manual(values = rev(cols), name = "Mean GMIS", labels = c("low (-1 sd)","high (+1 sd)")) +
  scale_fill_manual(values = rev(cols), name = "Mean GMIS", labels = c("low (-1 sd)","high (+1 sd)")) +
  # scale_y_continuous(breaks = seq(from = -0.4, to = 1.0, by = 0.2)) +
  ng1 + theme(legend.position = 'right')
LST_wNDVI_plot

ggsave(filename = "analysis/figures/manuscript-panels/figure-4/figure4C_HCNslope_wNDVIslope_sLST.pdf", 
       plot = LST_wNDVI_plot, device = "pdf", width = 8, height = 8, units = "in", dpi = 600, useDingbats = FALSE)



