# Script to create table with summary of HCN clines by city
# Uses OLS regression and not Robust regression
#
# Author: James S. Santangelo

# Create list with paths to dataframes
inpath <- "data/clean/popMeans_allCities_withEnviro/"
df_list <- create_df_list(inpath) %>% 
  map(., std_var_zero_one, var = 'GMIS_Mean')

# Get continent for each city
continents <- do.call(rbind, df_list) %>%
  group_by(city) %>%
  distinct(continent)

# Get stats from logistic regression by city using standardized distance as a predictor
log_reg_stats_distance <- df_list %>% 
  map_dfr(., logistic_regression_stats, 
          predictor = 'std_distance') %>% 
  rename_at(vars(-city),function(x) paste0(x,"_Dist"))

# Get change in log-odds of HCN from binomial regression
# These are extracted from the mixed model
logOdds_fromGlobalModel_forSupMat <- coef(glueClineModel_stdDist)$city %>%
  rownames_to_column(var = 'city') %>% 
  dplyr::select(city, std_distance) %>% 
  rename('betaLog_fromGlobalModel' = 'std_distance') %>% 
  mutate(betaLog_fromGlobalModel = round(betaLog_fromGlobalModel, 3))

# Get stats from logistic regression by city using standardized GMIS as a predictor
log_reg_stats_gmis <- df_list %>% 
  map_dfr(., logistic_regression_stats, 
          predictor = 'std_GMIS_Mean') %>% 
  rename_at(vars(-city),function(x) paste0(x,"_GMIS"))

# Get stats from logistic regression by city using standardized HII as a predictor
log_reg_stats_hii <- df_all_popMeans_stdHII %>% 
  filter(!(city == 'Elmira')) %>% 
  group_split(city) %>% 
  map_dfr(., logistic_regression_stats, 
          predictor = 'std_hii') %>% 
  rename_at(vars(-city),function(x) paste0(x,"_hii"))

# Merge logistic regression stats
linearClineTable_mod <- log_reg_stats_distance %>% 
  left_join(., log_reg_stats_gmis, by = 'city') %>% 
  left_join(., log_reg_stats_hii, by = 'city') %>% 
  left_join(., continents, by = 'city') %>% 
  left_join(., logOdds_fromGlobalModel_forSupMat, by = 'city')

# Write cline model summary to disk
outpath <- "analysis/tables/allCities_logisticReg_coefs.csv"
write_csv(linearClineTable_mod, outpath, col_names = TRUE)

