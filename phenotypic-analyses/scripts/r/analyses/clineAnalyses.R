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
inpath <- "data/clean/popMeans_allCities/"
popMeans_dfList <- create_df_list(inpath)
df_all_popMeans <- do.call(rbind, popMeans_dfList)

#########################################################
#### GLOBAL MODEL TESTING FOR THE EFFECT OF DISTANCE ####
#########################################################

# Model with standardized distance to get stats
glueClineModel_stdDist <- lmer(freqHCN ~ std_distance + continent + std_distance:continent + 
                         (1 + std_distance|city),
                       REML = FALSE, data = df_all_popMeans)
glueClineModel_stdDist_summary <- summary(glueClineModel_stdDist)
glueClineModel_stdDist_anova <- anova(glueClineModel_stdDist, type = 3)
glueClineModel_stdDist_r_squared <- MuMIn::r.squaredGLMM(glueClineModel_stdDist)

# Model with unstandardized distance to get average increase in HCN 
glueClineModel_unsStdDist <- lmer(freqHCN ~ distance + continent + distance:continent + 
                         (1 + std_distance|city),
                       REML = FALSE, data = df_all_popMeans)
# summary(glueClineModel)
glueClineModel_unsStdDist_distEffect <- coef(summary(glueClineModel_unsStdDist))["distance", "Estimate"]
avg_increaseHCN_perKM <- glueClineModel_unsStdDist_distEffect * 100


# Model with standardized distance and REML = TRUE to get stats for random effects
glueClineModel_randEff <- lmer(freqHCN ~ std_distance + continent + std_distance:continent + 
                         (1 + std_distance|city),
                       REML = TRUE, data = df_all_popMeans)
glueClineModel_randEffect_anova <- ranova(glueClineModel_randEff, reduce.terms = TRUE)
