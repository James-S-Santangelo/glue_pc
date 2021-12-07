# Script to look at effect of city metrics (e.g., Age, density, etc.) on HCN clines
#
# Author: James S. Santangelo
# Adapted from code by Marc Johnson and Hayden Fargo

###############
#### SETUP ####
###############

# Load data
city_stats <- read_csv('data/raw/city_data/City_characteristics.csv') %>% 
  rename('city' = 'City') %>% 
  dplyr::select(city, area, pop_size, density, no_cities, city_age, betaLog_Dist)

################
#### MODELS ####
################

# Model for effects of city stats on HCN using raw data
# Diagnostics suggest heterogeneity of variance and right skewed residuals
modraw <- lm(betaLog_Dist ~ area + pop_size + density + no_cities + city_age, data = city_stats)
# plot(modraw)
# hist(residuals(modraw))
summary(modraw)

# Log transform predictor variables
city_stats <- city_stats %>% 
  mutate(log_area = log(area),
         log_pop_size = log(pop_size),
         log_density = log(density),
         log_city_age = log(city_age + 1))

# Look at correlations of raw and log transformed variables
logVarsMat <- city_stats %>% dplyr::select(starts_with('log_')) %>% as.matrix()
rawVarsMat <- city_stats %>% dplyr::select(area, pop_size, density, no_cities, city_age) %>% as.matrix()
logVars_corr <- generate_correlation_df(logVarsMat)
rawVars_corr <- generate_correlation_df(rawVarsMat)

# Create the plots
pairs(logVarsMat, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
pairs(rawVarsMat, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

# Model with log-transformed predictors
modlog_all <- lm(betaLog_Dist ~ scale(log_area) + 
               scale(log_density) + 
               scale(no_cities) + 
               scale(log_city_age) +
               scale(log_pop_size), 
             data = city_stats)
# modlogden<-lm(hcn_slope~logdensity) #ran this because NaN comes out for logdensity on summary of full model
summary(modlog_all)
# plot(modlog_all)
# hist(residuals(modlog_all))

# Remove log_pop_size because colinear with area and age
# Diagnostics look good
modlog_popout<-lm(betaLog_Dist ~ scale(log_area) + 
                    scale(log_density) + 
                    scale(no_cities) + 
                    scale(log_city_age), 
                  data = city_stats)
# plot(modlog_popout)
# hist(residuals(modlog_popout))
summary(modlog_popout)

# Remove log_area because colinear with age
# Diagnostics look good
modlog_popout_areaout <- lm(betaLog_Dist ~ scale(log_density) + 
                    scale(no_cities) + 
                    scale(log_city_age), 
                  data = city_stats)
# plot(modlog_popout)
# hist(residuals(modlog_popout))
tableS7A_model <- summary(modlog_popout_areaout)

# Remove age because colinear with area 
# Diagnostics look good
modlog_popout_ageout <- lm(betaLog_Dist ~ scale(log_area) + 
                    scale(log_density) + 
                    scale(no_cities), 
                  data = city_stats)
# plot(modlog_popout)
# hist(residuals(modlog_popout))
tableS7B_model <- summary(modlog_popout_ageout)

#Filter out all city ages = 0
#remove logpop because colinear with area and age
#subset data to only include city_age>0
city_stats_ageAboveZero <- city_stats %>% 
  filter(city_age > 0)

# Model with denisty, no_cities, and age, excluding cities with estimated age = 0
modlog_ageAboveZero_areaout_popout <- lm(betaLog_Dist ~ log_density + 
                                           scale(no_cities) + 
                                           scale(log_city_age), 
                                         data = city_stats_ageAboveZero)
tableS7C_model <- summary(modlog_ageAboveZero_areaout_popout)

# Model with denisty, no_cities, and area, excluding cities with estimated age = 0
modlog_ageAboveZero_ageout_popout <- lm(betaLog_Dist ~ scale(log_density) + 
                                           scale(log_no_cities) + 
                                           scale(log_area), 
                                         data = city_stats_ageAboveZero)
summary(modlog_ageAboveZero_ageout_popout)
