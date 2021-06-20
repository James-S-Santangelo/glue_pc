# Script to look at effect of city metrics (e.g., Age, density, etc.) on HCN clines
#
# Author: James S. Santangelo
# Adapted from code by Marc Johnson and Hayden Fargo

###############
#### SETUP ####
###############

# Load data
city_stats <- read_csv('data/raw/city_data/City_characteristics.csv') %>% 
  dplyr::select(City, area, pop_size, density, no_cities, city_age, hcn_slope)

################
#### MODELS ####
################

# Model for effects of city stats on HCN using raw data
# Diagnostics suggest heterogeneity of variance and right skewed residuals
modraw <- lm(hcn_slope ~ area + pop_size + density + no_cities + city_age, data = city_stats)
# plot(modraw)
# hist(residuals(modraw))
summary(modraw)

# Log transform predictor variables
city_stats <- city_stats %>% 
  mutate(log_area = log(area),
         log_pop_size = log(pop_size),
         log_density = log(density),
         log_no_cities = log(no_cities),
         log_city_age = log(city_age + 1))

# Model with log-transformed predictors
modlog <- lm(hcn_slope ~ scale(log_area) + 
               scale(log_pop_size) + 
               scale(log_density) + 
               scale(log_no_cities) + 
               scale(log_city_age), 
             data = city_stats)
# modlogden<-lm(hcn_slope~logdensity) #ran this because NaN comes out for logdensity on summary of full model
summary(modlog)
# plot(modlog)
# hist(residuals(modlog))

# Remove log_pop_size because colinear with area and age
# Diagnostics look good
modlog_popout<-lm(hcn_slope ~ scale(log_area) + 
                    scale(log_density) + 
                    scale(log_no_cities) + 
                    scale(log_city_age), 
                  data = city_stats)
# plot(modlog_popout)
# hist(residuals(modlog_popout))
finalModel_cityStats <- summary(modlog_popout)
