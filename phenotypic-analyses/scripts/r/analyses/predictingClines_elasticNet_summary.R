# Script to summarize elastic net model coefficients and perform PC regression
#
# Author: James S. Santangelo

############################################
#### SUMMARIZE ELASTIC NET COEFFICIENTS ####
############################################

# Load in Elastic Net coefficients and results
elasticNet_obs_coefs <- read_csv('analysis/supplementary-tables/elasticNet_obs_coefs.csv')
elasticNet_obs_results <- read_csv('analysis/supplementary-tables/elasticNet_obs_result.csv')

# Extract nonzero coefficients
elasticNet_obs_coefs_nonZeroOnly <- elasticNet_obs_coefs %>% 
  dplyr::select(-index) %>% 
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0)
elasticNet_obs_coefs_names <- colnames(elasticNet_obs_coefs)

# Summarize coefficients
elasticNet_coefSummary <- elasticNet_obs_coefs_nonZeroOnly %>% 
  pivot_longer(names(.), names_to = 'predictor', values_to = 'coef') %>% 
  group_by(predictor) %>%
  mutate(is_non_zero = ifelse(coef == 0, 0, 1)) %>% 
  summarise(mean = mean(coef),
            non_zero = sum(is_non_zero)) %>% 
  arrange(desc(non_zero))

write_csv(elasticNet_coefSummary, 'analysis/supplementary-tables/elasticNet_coefSummary.csv')