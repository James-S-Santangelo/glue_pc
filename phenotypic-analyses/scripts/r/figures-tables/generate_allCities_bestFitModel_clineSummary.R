# Script to create table with summary of HCN clines by city
# Uses OLS regression and not Robust regression
#
# Author: James S. Santangelo

# Create list with paths to dataframes
inpath <- "data/clean/popMeans_allCities_withEnviro/"
df_list <- create_df_list(inpath)

# Get continent for each city
continents <- do.call(rbind, df_list) %>%
  group_by(city) %>%
  distinct(continent)

# Add continent to summary of linear models for clines
# Model stats are from best fit model ('linear' or 'quadratic')
linearClineTable_mod <- clineResults(df_list) %>%
  left_join(., continents %>% dplyr::select(continent), by = 'city')

# Get stats from logistic regression by city
log_reg_stats <- df_list %>% 
  map_dfr(., logistic_regression_stats)

# Merge logistic regression stats
linearClineTable_mod <- linearClineTable_mod %>% 
  left_join(., log_reg_stats, by = 'city')

# Write cline model summary to disk
outpath <- "analysis/supplementary-tables/allCities_bestFitModel_clineSummary.csv"
write_csv(linearClineTable_mod, outpath, col_names = TRUE)

