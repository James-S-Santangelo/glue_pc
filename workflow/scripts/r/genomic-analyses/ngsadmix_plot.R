# Script to generate admixture plots


# Load sample order
habitat_info <- read_delim('resources/sequencedPlants_phenotypesHabitat.txt', 
                           delim = '\t') %>% 
  dplyr::select(Sample, Habitat, Population, Plant)
samples <- read_table('data/angsd_sample_order.txt', col_names = FALSE) %>% 
  rename('Sample' = 'X1') %>%
  left_join(., habitat_info, by = 'Sample')

# K = 2
k2 <- read_delim('data/ngsadmix/K2/allSamples_ngsadmix_4fold_maf0.05_K2_seed1.qopt',
                 col_names = FALSE, delim = ' ') %>% 
  bind_cols(., samples) %>% 
  dplyr::select(-X3) %>% 
  pivot_longer(X1:X2, values_to = 'Probs') %>% 
  mutate(Probs = round(Probs, 5)) %>% 
  arrange(Habitat, Population, Plant)

cols <- wes_palette("Darjeeling1", n = 2, type = 'discrete')
k2plot <-
  ggplot(k2, aes(factor(Sample), Probs, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Habitat), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_manual(values = cols) + 
  theme(
    legend.position = 'none',
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(17,1,1,1),"cm")
  ) 
k2plot

# K = 3
k3 <- read_delim('data/ngsadmix/K3/allSamples_ngsadmix_4fold_maf0.05_K3_seed1.qopt',
                 col_names = FALSE, delim = ' ') %>% 
  bind_cols(., samples) %>% 
  dplyr::select(-X4) %>% 
  pivot_longer(X1:X3, values_to = 'Probs') %>% 
  mutate(Probs = round(Probs, 5)) %>% 
  arrange(Habitat, Population, Plant)

cols <- wes_palette("Darjeeling1", n = 3, type = 'discrete')
k3plot <-
  ggplot(k3, aes(factor(Sample), Probs, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Habitat), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_manual(values = cols) + 
  theme(
    legend.position = 'none',
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(17,1,1,1),"cm")
  ) 
k3plot

# K = 4
k4 <- read_delim('data/ngsadmix/K4/allSamples_ngsadmix_4fold_maf0.05_K4_seed1.qopt',
                 col_names = FALSE, delim = ' ') %>% 
  bind_cols(., samples) %>% 
  dplyr::select(-X5) %>% 
  pivot_longer(X1:X4, values_to = 'Probs') %>% 
  mutate(Probs = round(Probs, 5)) %>% 
  arrange(Habitat, Population, Plant)

cols <- wes_palette("Darjeeling1", n = 4, type = 'discrete')
k4plot <-
  ggplot(k4, aes(factor(Sample), Probs, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Habitat), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_manual(values = cols) + 
  theme(
    legend.position = 'none',
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(17,1,1,1),"cm")
  ) 
k4plot

# K = 5
k5 <- read_delim('data/ngsadmix/K5/allSamples_ngsadmix_4fold_maf0.05_K5_seed1.qopt',
                 col_names = FALSE, delim = ' ') %>% 
  bind_cols(., samples) %>% 
  dplyr::select(-X6) %>% 
  pivot_longer(X1:X5, values_to = 'Probs') %>% 
  mutate(Probs = round(Probs, 5))

cols <- wes_palette("Darjeeling1", n = 5, type = 'discrete')
k5plot <-
  ggplot(k5, aes(factor(Sample), Probs, fill = factor(name))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(Habitat), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_fill_manual(values = cols) + 
  theme(
    legend.position = 'none',
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) 
k5plot
