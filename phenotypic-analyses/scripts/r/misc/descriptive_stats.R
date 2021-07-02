# Script to calculate some descriptive statistics for the manuscript
#
# Author: James S. Santangelo

# Mean number of plants per population with standard errors
mean_plants_per_pop <- allPopMeans %>% 
  mutate(total_plants = if_else(is.na(total_plants), 20, total_plants)) %>% 
  summarise(mean = mean(total_plants),
            total_pops = n(),
            se = mean / sqrt(total_pops),
            min = min(total_plants), 
            max = max(total_plants))

# Mean number of populations per city
mean_populations <- allPopMeans %>% 
  group_by(city) %>% 
  summarise(num_pops = n()) %>% 
  ungroup() %>% 
  summarise(mean_pops = mean(num_pops),
            num_cities = n(),
            se = mean_pops / sqrt(num_cities))

# Total number of plants
num_plants <- allPopMeans %>% 
  mutate(total_plants = if_else(is.na(total_plants), 20, total_plants)) %>% 
  summarise(num_plants = sum(total_plants)) %>% pull()

# Total number of populations
# Each row in global table is a population
num_populations <- allPopMeans %>% nrow()

# Total number of cities
num_cities <- final_table %>% nrow()

# Percent significant clines from indipendent logistic regressions and by Distance or GMIS as predictor
percent_sig_clines_logReg <- linearClineTable_mod %>% 
  dplyr::select(city, pvalLog_Dist, pvalLog_GMIS) %>% 
  pivot_longer(cols = c('pvalLog_Dist', 'pvalLog_GMIS'), names_to = 'group', values_to = 'pvalLog') %>% 
  mutate(sig = case_when(pvalLog < 0.05 ~ 'Yes',
                         TRUE ~ 'No')) %>% 
  group_by(sig, group) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)

# Percent significant by direction and by Distance or GMIS as predictor
percent_sig_clines_logReg_byDirection <- linearClineTable_mod %>% 
  dplyr::select(city, pvalLog_Dist, pvalLog_GMIS, betaLog_Dist, betaLog_GMIS) %>% 
  pivot_longer(cols = c('pvalLog_Dist', 'pvalLog_GMIS'), names_to = 'group_sig', values_to = 'pvalLog') %>% 
  pivot_longer(cols = c('betaLog_Dist', 'betaLog_GMIS'), names_to = 'group_dir', values_to = 'betaLog') %>% 
  mutate(sig = case_when(pvalLog < 0.05 ~ 'Yes',
                         TRUE ~ 'No'),
         direction = case_when(betaLog < 0 ~ 'Neg',
                               TRUE ~ 'Pos')) %>% 
  group_by(group_sig, group_dir, direction, sig) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100) %>% 
  filter((group_sig  == 'pvalLog_Dist' & group_dir == 'betaLog_Dist') |
          (group_sig  == 'pvalLog_GMIS' & group_dir == 'betaLog_GMIS'))

