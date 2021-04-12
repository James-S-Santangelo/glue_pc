# Script to perform environmental analyses for the first question of the paper:
#   Does urbanization lead to convergent environmental change in cities throughout the world?
#
# Authors: James S. Santangelo and Pedro Peres-Neto

###############
#### SETUP ####
###############

# Create dataframe with population-mean HCN and environmental variables for every city
# Needs to be created using `read.csv()` for Pedro's code to work...
inpath <- "data/clean/popMeans_allCities_withEnviro/"
csv.files <- list.files(path = inpath, pattern="*.csv")
df_all_popMeans <- c()
for (i in 1:length(csv.files)){
  data <- read.csv(paste0(inpath, csv.files[i])) %>% dplyr::select(city, 
                                                                   std_distance, 
                                                                   freqHCN, 
                                                                   population_latitude,
                                                                   population_longitude,
                                                                   matches("*Mean$"))
  df_all_popMeans <- rbind(df_all_popMeans, data)
}

#######################################################
#### CONVERGENCE IN MEAN MULTIVARIATE ENVIRONMENTS ####
#######################################################

# `generate.pred.values` performs a robuts regression with mean environmental values as the response variable and standardised distance
# to the urban core as a predictor. Done for each variable and city separately. 
enviroRLM_results <- generate.pred.values(all.data = df_all_popMeans, permute=FALSE)

# Urban/rural predicted values of environmental variables
Predicted_vals<- enviroRLM_results$Predicted.Values

# Extreme values of urban and rural environmental variables from predicted dataset
Extreme_vals <- pick.extreme.values(Predicted.Values = Predicted_vals, Original.Values = df_all_popMeans, number.extreme.sites = 1)
UrbanPredExtremes <- Extreme_vals$UrbanPredExtremes %>% as.data.frame() %>% dplyr::select(-freqHCN) %>% as.matrix()
RuralPredExtremes <- Extreme_vals$RuralPredExtremes %>% as.data.frame() %>% dplyr::select(-freqHCN) %>% as.matrix()
Extreme_vals$city.names
# Vector of city names
city_names <- unique(Extreme_vals$city.names)
n_cities <- length(unique(city_names))
habitat <- c(rep("Urban", n_cities), rep("Rural", n_cities))

# PCA of urban and rural extreme environmental variables
# RDA from vegan package perform PCA if no predictor is provided
enviroPCA <- rda(rbind(UrbanPredExtremes, RuralPredExtremes), 
                 scale = TRUE, na.action = "na.omit")
enviroPCA_summary <- summary(enviroPCA)

# Permutation test for difference in mean environmental variables between urban and rural habitats
# Can only test for interaction if number.extreme.sites >= 2
enviroPCA_permute <- permutation.tests(all.data = df_all_popMeans, n.perm = n_perm, number.extreme.sites = 2)
# 
# ########################################################
# #### ENVIRONMENTAL VARIANCE IN URBAN/RURAL HABITATS ####
# ########################################################
# 
# Urban-rural variance in environmental variables between cities
res.dist <- mult.dispersion(all.data = df_all_popMeans, number.extreme.sites = 5)
enviroVariance <- multi.disp.analysis(res.dist = res.dist, plot.disp = FALSE)

# BoxM test for homogeneity between two covariance matrices (i.e., urban and rural predicted environmental variable matrices.)
enviroVariance_boxM <- enviroVariance$boxM

# Multivariate levene's test for homogeneity of urban and rural environmental variances
enviroVariance_leveneLM <- enviroVariance$levene_lm
enviroVariance_leveneLM_anova <- Anova(enviroVariance$levene_lm)
enviroVariance_leveneCAN <- enviroVariance$levene_can

