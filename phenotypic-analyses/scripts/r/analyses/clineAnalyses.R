# Main script for analyzing clines in HCN across cities for GLUE
#
# Author: James S. Santangelo
# Date: Saturday Jan. 12, 2019

###############
#### SETUP ####
###############

# Change default contrasts to enable type III SS
options(contrasts = c("contr.sum", "contr.poly"))
# options(contrasts = c("contr.treatment", "contr.poly")) # Default

# Dataframe with population-Mean HCN frequencies for all cities.
# Includes continent and country columns.
inpath <- "data/clean/popMeans_allCities_withEnviro/"
popMeans_dfList <- create_df_list(inpath)
df_all_popMeans <- do.call(rbind, popMeans_dfList) %>% 
  filter(!(is.na(std_distance)))

#########################################################
#### GLOBAL MODEL TESTING FOR THE EFFECT OF DISTANCE ####
#########################################################

## MODEL WITH STANDARDIZED DISTANCE ##

# Generalized linear model with population-mean HCN frequencies as response and standardized distance
# to the urban core as predictor (+ continent and distance:continent interactions). Binomial error
# distribution with per-population sample sizes provided as weights.

# Random-slope model allowing effect of distance to vary in each city
# Estimates correlation between random intercept and slope
glueClineModel_stdDist <- glmer(freqHCN ~ std_distance + continent + std_distance:continent + 
                                  (1 + std_distance |city),
                                data = df_all_popMeans,
                                weights = total_plants,
                                family = 'binomial',
                                control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
glueClineModel_stdDist_summary <- summary(glueClineModel_stdDist)
glueClineModel_stdDist_anova <- Anova(glueClineModel_stdDist, type = 3)
glueClineModel_stdDist_r_squared <- MuMIn::r.squaredGLMM(glueClineModel_stdDist)

# Model with no random effects
glueClineModel_stdDist_noRE <- glm(freqHCN ~ std_distance + continent + std_distance:continent,
                              data = df_all_popMeans,
                              weights = total_plants,
                              family = 'binomial')

# LRT test to see if estimating correlation improves fit
glueClineModel_stdDist_noCor <- update(glueClineModel_stdDist, . ~ . -(1 + std_distance|city) + (1 | city) + (0 + std_distance | city))
glueClineModel_stdDist_corRE_test <- anova(glueClineModel_stdDist_noCor, glueClineModel_stdDist)

# LRT test to see if random slope model better fit than intercept-only model
glueClineModel_stdDist_intOnly <- update(glueClineModel_stdDist_noCor, . ~ . -(0 + std_distance | city))
glueClineModel_stdDist_slopeRE_test <- anova(glueClineModel_stdDist_intOnly, glueClineModel_stdDist_noCor)

# LRT to see if estimating random intercepts improves fit
glueClineModel_stdDist_randInt_test <- anova(glueClineModel_stdDist_intOnly, glueClineModel_stdDist_noRE)

# Predicted proportion of HCN+ across transect from mixed-model
propHCN_predicted <- ggeffects::ggeffect(glueClineModel_stdDist, 
                                         terms = c('std_distance[all]'), 
                                         type = 're')

dist0_pred <- propHCN_predicted %>% filter(x == 0) %>% pull(predicted)
dist1_pred <- propHCN_predicted %>% filter(x == 1) %>% pull(predicted)
pred_percent_change <- round((dist1_pred - dist0_pred) / dist0_pred, 2) * 100

###################################################
#### SUPPLEMENTARY ANALYSES: GMIS AS PREDICTOR ####
###################################################

# Re-run same analysis as above but using standardized GMIS as a predictor instead of 
# standardized distance

# Create dataframe with standardized GMIS as column
df_all_popMeans_stdGMIS <- df_all_popMeans %>% 
  group_split(city) %>% 
  map_dfr(., std_var_zero_one, var = 'GMIS_Mean')

# Random-slope model allowing effect of GMIS to vary in each city
# Estimates correlation between random intercept and slope
glueClineModel_gmis <- glmer(freqHCN ~ std_GMIS_Mean + continent + std_GMIS_Mean:continent + 
                               (1 + std_GMIS_Mean |city),
                             data = df_all_popMeans_stdGMIS,
                             weights = total_plants,
                             family = 'binomial',
                             control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
glueClineModel_stdGmis_summary <- summary(glueClineModel_gmis)
glueClineModel_stdGmis_anova <- Anova(glueClineModel_gmis, type = 3)
glueClineModel_stdGmis_r_squared <- MuMIn::r.squaredGLMM(glueClineModel_gmis)

# Model with no random effects
glueClineModel_gmis_noRE <- glm(freqHCN ~ std_GMIS_Mean + continent + std_GMIS_Mean:continent,
                             data = df_all_popMeans_stdGMIS,
                             weights = total_plants,
                             family = 'binomial')

# LRT test to see if estimating correlation improves fit
glueClineModel_stdGmis_noCor <- update(glueClineModel_gmis, . ~ . -(1 + std_GMIS_Mean|city) + (1 | city) + (0 + std_GMIS_Mean | city))
glueClineModel_stdGmis_corRE_test <- anova(glueClineModel_stdGmis_noCor, glueClineModel_gmis)

# LRT test to see if random slope model better fit than intercept-only model
glueClineModel_stdGmis_intOnly <- update(glueClineModel_stdGmis_noCor, . ~ . -(0 + std_GMIS_Mean | city))
glueClineModel_stdGmis_slopeRE_test <- anova(glueClineModel_stdGmis_intOnly, glueClineModel_stdGmis_noCor)

# LRT to see if estimating random intercepts improves fit
glueClineModel_stdGmis_randInt_test <- anova(glueClineModel_stdGmis_intOnly, glueClineModel_gmis_noRE)

##################################################
#### SUPPLEMENTARY ANALYSES: HII AS PREDICTOR ####
##################################################

# Re-run same analysis as above but using standardized Human Influence Index as a predictor instead of 
# standardized distance

# Create dataframe with HII values for all cities
hii_inpath <- 'data/raw/environmental_data/hii/'
hii_df <- create_df_list(hii_inpath) %>%
  do.call(rbind, .) %>% 
  dplyr::select(city, population, hii) %>% 
  mutate(city = fct_recode(city, 'Newhaven' = 'New_Haven'))

# Join HII values to population-mean dataframes and standardize
df_all_popMeans_stdHII <- df_all_popMeans %>% 
  left_join(., hii_df, by = c("city", "population")) %>% 
  group_split(city) %>% 
  map(., std_var_zero_one, var = 'hii') %>% 
  do.call(rbind, .)

# Random-slope model allowing effect of HII to vary in each city
# Estimates correlation between random intercept and slope
glueClineModel_hii <- glmer(freqHCN ~ std_hii + continent + std_hii:continent + 
                               (1 + std_hii |city),
                             data = df_all_popMeans_stdHII,
                             weights = total_plants,
                             family = 'binomial',
                             control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
glueClineModel_stdhii_summary <- summary(glueClineModel_hii)
glueClineModel_stdhii_anova <- Anova(glueClineModel_hii, type = 3)
glueClineModel_stdhii_r_squared <- MuMIn::r.squaredGLMM(glueClineModel_hii)

# Model with no random effects
glueClineModel_hii_noRE <- glm(freqHCN ~ std_hii + continent + std_hii:continent,
                            data = df_all_popMeans_stdHII,
                            weights = total_plants,
                            family = 'binomial')

# LRT test to see if estimating correlation improves fit
glueClineModel_stdhii_noCor <- update(glueClineModel_hii, . ~ . -(1 + std_hii|city) + (1 | city) + (0 + std_hii | city))
glueClineModel_stdhii_corRE_test <- anova(glueClineModel_stdhii_noCor, glueClineModel_hii)

# LRT test to see if random slope model better fit than intercept-only model
glueClineModel_stdhii_intOnly <- update(glueClineModel_stdhii_noCor, . ~ . -(0 + std_hii | city))
glueClineModel_stdhii_slopeRE_test <- anova(glueClineModel_stdhii_intOnly, glueClineModel_stdhii_noCor)

# LRT to see if estimating random intercepts improves fit
glueClineModel_stdhii_randInt_test <- anova(glueClineModel_stdhii_intOnly, glueClineModel_hii_noRE)
