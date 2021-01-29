# Script to examine effects of BAQ and GL model on SFS

# Load required packages
library(tidyverse)
library(wesanderson)

###################
#### FUNCTIONS ####
###################

load_wide_sfs <- function(path){
  
  # Get name of folder with parameter combinations
  dirname <- dirname(path)
  
  # Read in SFS
  full_path <- paste0(inpath, '/', path)
  sfs <- read_table(full_path, col_names = FALSE) %>% 
    t() %>% 
    as.data.frame() %>% 
    rename('num_sites' = 'V1') %>% 
    mutate(bin = 1:n() - 1,
           id = dirname,
           prop_sites = num_sites / sum(num_sites),
           id = dirname)%>% 
    filter(bin <= num_samples)
  return(sfs)
}

#####################
#### COMPARE BAQ ####
#####################

# Load in all data
inpath <- 'data/test_params/'
num_samples <- 120
sfs_df <- list.files(inpath, recursive = TRUE) %>% 
  map_dfr(., load_wide_sfs) %>% 
  separate(id, into = c('gl', 'baq'), sep = '_')

# Samtools SFS with 3 different BAQ algorithms
# Invariant sites not plotted
# Only show up to 25-tons
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
sfs_df %>% 
  filter(bin != 0 & bin <= 25)  %>%
  filter(gl == 'GL1') %>%
  ggplot(., aes(x = bin, y = prop_sites, fill = baq)) + 
    geom_bar(stat ='identity', color = 'black',  width=.70, position = "dodge") + 
  ylab('Proportion of sites') + xlab('Allele frequency') +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(1, 40, 3)) +
  scale_y_continuous(breaks = seq(0, 0.035, 0.005)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))
  

# GATK SFS with 3 different BAQ algorithms
# Invariant sites not plotted
# Only show up to 25-tons
sfs_df %>% 
  filter(bin != 0 & bin <= 25)  %>%
  filter(gl == 'GL2') %>%
  ggplot(., aes(x = bin, y = prop_sites, fill = baq)) + 
  geom_bar(stat ='identity', color = 'black',  width=.70, position = "dodge") + 
  ylab('Proportion of sites') + xlab('Allele frequency') +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = seq(1, 40, 3)) +
  scale_y_continuous(breaks = seq(0, 0.045, 0.005)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

