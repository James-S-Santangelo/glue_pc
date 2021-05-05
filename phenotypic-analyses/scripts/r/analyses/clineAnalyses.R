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

# LRT test to see if estimating correlation improves fit
glueClineModel_stdDist_noCor <- update(glueClineModel_stdDist, . ~ . -(1 + std_distance|city) + (1 | city) + (0 + std_distance | city))
glueClineModel_stdDist_corRE_test <- anova(glueClineModel_stdDist_noCor, glueClineModel_stdDist)

# LRT test to see if random slope model better fit than intercept-only model
glueClineModel_stdDist_intOnly <- update(glueClineModel_stdDist_noCor, . ~ . -(0 + std_distance | city))
glueClineModel_stdDist_slopeRE_test <- anova(glueClineModel_stdDist_intOnly, glueClineModel_stdDist_noCor)

# LRT test between model with random slope and correlation, and intercept only
glueClineModel_stdDist_slopeREwithCor_test <- anova(glueClineModel_stdDist_intOnly, glueClineModel_stdDist)

# Predicted proportion of HCN+ across transect from mixed-model
propHCN_predicted <- ggeffects::ggeffect(glueClineModel_stdDist, 
                                         terms = c('std_distance[all]'), 
                                         type = 're')

