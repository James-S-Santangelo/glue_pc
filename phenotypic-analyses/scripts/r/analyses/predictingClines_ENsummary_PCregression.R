# Script to summarize elastic net model coefficients and perform PC regression
#
# Author: James S. Santangelo

# Load in Elastic Net coefficients and results
num_reps <- 100
elasticNet_obs_coefs <- read_csv('analysis/supplementary-tables/elasticNet_obs_coefs.csv')
elasticNet_obs_results <- read_csv('analysis/supplementary-tables/elasticNet_obs_result.csv')

elasticNet_obs_coefs_nonZeroOnly <- elasticNet_obs_coefs %>% 
  dplyr::select(-index) %>% 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0)
elasticNet_obs_coefs_names <- colnames(elasticNet_obs_coefs)

elasticNet_coefSummary <- elasticNet_obs_coefs_nonZeroOnly %>% 
  pivot_longer(names(.), names_to = 'predictor', values_to = 'coef') %>% 
  group_by(predictor) %>%
  mutate(is_non_zero = ifelse(coef == 0, 0, 1)) %>% 
  summarise(mean = mean(coef),
            non_zero = sum(is_non_zero)) %>% 
  arrange(desc(non_zero))

write_csv(elasticNet_coefSummary, 'analysis/supplementary-tables/elasticNet_coefSummary.csv')


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



