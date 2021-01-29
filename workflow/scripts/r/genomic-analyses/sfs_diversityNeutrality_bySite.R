# Script to plot SFS for all Sites and 4 fold sites and calculate diversity/neutrality stats

# Load packages
library(tidyverse)

###################
#### FUNCTIONS ####
###################

load_long_sfs <- function(path){
  
  # Get name of folder with parameter combinations
  dirname <- dirname(path)
  
  # Read in SFS
  full_path <- paste0(inpath, '/', path)
  sfs <- read_delim(full_path, delim= '\t', col_names = FALSE) %>% 
    as.data.frame() %>% 
    rename('bin' = 'X1',
           'num_sites' = 'X2') %>% 
    filter(bin <= num_samples) %>% 
    mutate(prop_sites = num_sites / sum(num_sites),
           site = dirname)
  return(sfs)
}

load_div_neut_df <- function(path){
  
  # Get name of folder with parameter combinations
  dirname <- dirname(path)
  
  # Read in stats
  full_path <- paste0(inpath, '/', path)
  stats <- read_delim(full_path, delim= '\t', col_names = TRUE) %>% 
    mutate(site = dirname,
           tp_scaled = tP / nSites,
           tw_scaled = tW / nSites)
  return(stats)
  
}


##################
#### PLOT SFS ####
##################

inpath <- 'data/sfs/'
num_samples <- 120
sfs_df <- list.files(inpath, recursive = TRUE) %>% 
  map_dfr(., load_long_sfs)

# SFS with for all sites, 0fold, and 4fold sites
# Invariant sites not plotted
# Only show up to 25-tons
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
sfs_df %>% 
  filter(bin != 0 & bin <= 25)  %>%
  ggplot(., aes(x = bin, y = prop_sites, fill = site)) + 
  geom_bar(stat ='identity', color = 'black',  width=.70, position = "dodge") + 
  ylab('Proportion of sites') + xlab('Allele frequency') +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(1, 40, 3)) +
  scale_y_continuous(breaks = seq(0, 0.035, 0.005)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

##################################
#### DIVERSITY AND NEUTRALITY ####
##################################

inpath <- 'data/summary_stats/thetas/'
stats_df <- list.files(inpath, recursive = TRUE) %>% 
  map_dfr(., load_div_neut_df)

summary_stats <- stats_df %>% 
  group_by(site) %>% 
  summarise(total_sites = sum(nSites),
            mean_tp = mean(tp_scaled),
            mean_tw = mean(tw_scaled), 
            mean_td = mean(Tajima))

