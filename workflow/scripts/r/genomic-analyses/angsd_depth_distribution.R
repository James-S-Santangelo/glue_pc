# Script to look at depth distribution across chromosome 1

# Load packages
library(tidyverse)
library(wesanderson)

# Load depth dataset
num_samples = 120
set.seed(42)
depth_df_chrom1 <- read_delim('data/depth/CM019101.1/CM019101.1_sum_site_depths.txt', delim = '\t', 
                              col_names = FALSE) %>% 
  rename('total_depth' = 'X1') %>% 
  mutate(dp_per_ind = total_depth / num_samples) 

depth_df_chrom1_forPlot <- depth_df_chrom1 %>% 
  sample_frac(0.01)

ggplot(depth_df_chrom1_forPlot, aes(x = total_depth)) +
  geom_histogram(bins = 200, color = 'black', fill = 'white') + 
  xlab('Total depth') + ylab('Number of sites') +
  theme_classic() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

