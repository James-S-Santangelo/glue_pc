# Script to perform PCA on covariance in allele frequencies at 4-fold sites

# Load packages
library(tidyverse)
library(wesanderson)

# Load covariance matrix
cov_mat <- read_delim('data/pcangsd/allSamples_4fold_maf0.05_pcangsd.cov', 
                      col_names = FALSE, delim = ' ') %>% 
  as.matrix()

# Load sample order
habitat_info <- read_delim('resources/sequencedPlants_phenotypesHabitat.txt', 
                           delim = '\t') %>% 
  dplyr::select(Sample, Habitat, Population, Plant)
samples <- read_table('data/angsd_sample_order.txt', col_names = FALSE) %>% 
  rename('Sample' = 'X1') %>%
  left_join(., habitat_info, by = 'Sample')
  
# Extract eigenvectors and create df
eigenvectors <- eigen(cov_mat)
eigen_df <- eigenvectors$vectors %>% 
  as.data.frame() %>% 
  dplyr::select(V1, V2) %>% 
  rename('PC1' = 'V1',
         'PC2' = 'V2') %>% 
  bind_cols(., habitat_info)

# Plot
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
ggplot(eigen_df, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  xlab('PC1 (50%)') + ylab('PC2 (1.4%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))

#####################################
#### EXCLUDING WEIRD INDIVIDUALS ####
#####################################

problem_samples <- eigen_df %>% 
  filter(PC1 < -0.1 | PC2 < -0.3) %>% 
  pull(Sample)

problem_samples_indices <- samples %>% 
  mutate(index = 1:n()) %>% 
  filter(Sample %in% problem_samples) %>% 
  pull(index)

cov_mat_reduced <- cov_mat[-problem_samples_indices, -problem_samples_indices]
samples_good <- samples %>% 
  filter(!(Sample %in% problem_samples))


# Extract eigenvectors and create df
summary(princomp(cov_mat_reduced))
eigenvectors_reduced <- eigen(cov_mat_reduced)
eigen_df_reduced <- eigenvectors_reduced$vectors %>% 
  as.data.frame() %>% 
  dplyr::select(V1, V2) %>% 
  rename('PC1' = 'V1',
         'PC2' = 'V2') %>% 
  bind_cols(., samples_good)

# Plot
cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
ggplot(eigen_df_reduced, aes(x = PC1, y = PC2, color = Habitat)) +
  geom_point(size = 2.5) + 
  scale_color_manual(values = cols) + 
  theme_classic() + 
  xlab('PC1 (2.6%)') + ylab('PC2 (2.2%)') +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 15))
