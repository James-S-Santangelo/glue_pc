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

# Percent significant clines
percent_sig_clines_linOnly <- final_table %>% 
  group_by(sigRLM) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)

# How does this change with quadratic clines?
percent_sig_clines_withQuad <- linearClineTable_mod %>% 
  mutate(sig = case_when(modelOrderBestFit == 'linear' & pvalLin < 0.05 ~ 'Yes',
                         modelOrderBestFit == 'quadratic' & pvalLin < 0.05 & pvalQuad < 0.05 ~ 'Yes',
                         TRUE ~ 'No')) %>% 
  group_by(sig) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)

# Percent clines by direction. Linear RLM only
percent_sig_clines_byDirection <- final_table %>% 
  mutate(direction = ifelse(betaRLM_freqHCN < 0, 'Negative', 'Positive')) %>% 
  group_by(direction, sigRLM) %>% 
  summarise(count = n(),
            percent = (count / num_cities) * 100)

