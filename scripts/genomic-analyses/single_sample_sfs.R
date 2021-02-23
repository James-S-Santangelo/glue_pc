# Script to look at single sample SFS for all 500 samples to look for oddities

# Load required packages
library(tidyverse)

###################
#### FUNCTIONS ####
###################

load_wide_sfs <- function(path){
  
  # Get name of folder with parameter combinations
  dirname <- dirname(path)
  
  # Read in SFS
  full_path <- paste0(inpath, '/', path)
  sfs <- read_table(full_path, col_names = FALSE) %>% 
    as.data.frame() %>% 
    rename('Invar' = 'X1',
           'Var' = 'X2') %>% 
    dplyr::select(-X3) %>% 
    mutate(Sample = dirname,
           total_sites = Invar + Var,
           prop_var = Var / total_sites)

  return(sfs)
}


# Load in all data
inpath <- 'data/single_sample_sfs/'
sfs_df <- list.files(inpath, recursive = TRUE) %>% 
  map_dfr(., load_wide_sfs)

ggplot(data = sfs_df, aes(x = prop_var)) + 
  geom_histogram(bins = 100, color = 'black', fill = 'white') + 
  geom_vline(xintercept = 0.00155, linetype = 'dashed', color = 'red') +
  geom_vline(xintercept = 0.0103, linetype = 'dashed', color = 'red') +
  theme_classic() 

ggplot(data = sfs_df, aes(x = total_sites)) + 
  geom_histogram(bins = 100, color = 'black', fill = 'white') + 
  theme_classic()

sfs_df_red <- sfs_df %>% 
  filter(prop_var < 0.00155)

sfs_df_red <- sfs_df %>% 
  filter(prop_var > 0.0103)

#### MULTIQC DATA ####

multiqc <- read_delim('~/Sync/multiqc/multiqc_data/multiqc_qualimap_bamqc_genome_results_qualimap_bamqc.txt', 
                      delim='\t') 

multi_qc_split <- multiqc %>% 
  mutate(Sample = str_extract(Sample, pattern = '\\w+_\\d+_\\d+$')) %>% 
  dplyr::select(Sample, mean_coverage, percentage_aligned, general_error_rate) %>% 
  left_join(., sfs_df, by = 'Sample')

multi_qc_split %>% 
  ggplot(., aes(x = mean_coverage)) +
  geom_histogram(bins = 100, color = 'black', fill = 'white') + 
  geom_vline(xintercept = 0.35, color = 'red') +
  theme_classic()

ggplot(data = multi_qc_split, aes(x = mean_coverage, y = percentage_aligned, color = prop_var)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split, aes(x = mean_coverage, y = prop_var, color = percentage_aligned)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split, aes(x = percentage_aligned, y = prop_var, color = mean_coverage)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split, aes(x = general_error_rate, y = prop_var)) + 
  geom_point() +
  theme_classic()

multi_qc_split_red <- multi_qc_split %>% 
  filter(mean_coverage >= 0.31 & general_error_rate < 0.04)

bad_samples <- multi_qc_split %>% 
  filter(mean_coverage < 0.31 | general_error_rate > 0.04) %>% 
  left_join(dat, by = 'Sample') %>% 
  group_by(city, site) %>% 
  summarise(count = n()) %>% 
  as.data.frame()

dat <- read_csv('data-clean/low-1/low1_libraryConcentrations.csv') %>% 
  select(city, site, plantID) %>% 
  rename('Sample' = 'plantID')

write_csv(bad_samples, '~/Desktop/glue-low1_badSamples.csv')
  
multi_qc_split_red %>% 
  ggplot(., aes(x = mean_coverage)) +
  geom_histogram(bins = 100, color = 'black', fill = 'white') + 
  geom_vline(xintercept = 0.35, color = 'red') +
  theme_classic()

ggplot(data = multi_qc_split_red, aes(x = mean_coverage, y = percentage_aligned, color = prop_var)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split_red, aes(x = mean_coverage, y = prop_var, color = percentage_aligned)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split_red, aes(x = percentage_aligned, y = prop_var, color = mean_coverage)) + 
  geom_point() +
  theme_classic()

ggplot(data = multi_qc_split_red, aes(x = general_error_rate, y = prop_var)) + 
  geom_point() +
  theme_classic()
