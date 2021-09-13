# Script to create sample sheet for plants sequenced as part of DEEP3, LANE2

# Load in data and add native vs. introduced range
isdup <- function (x) duplicated (x)
deep3_low2_lane3_sampleSheet <- read_csv('data/clean/deep3/deep3_low2_lane3_dilutions.csv') %>% 

  dplyr::select(continent, city, pop, individual, lib_split, site, plantID, lane, i5_primer_index_seq, i7_primer_index_seq) %>% 
  rename('sample' = 'plantID') %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                           city == 'Tehran' ~ 'Native',
                           TRUE ~ 'Introduced')) %>% 
  
  # Add '_dup_' suffix to duplicated samples
  mutate(is_dup = isdup(sample),
         sample = ifelse(is_dup, paste0(sample, '_dup'), sample)) %>% 
  dplyr::select(-is_dup)
  

outpath <- 'resources/'
print(sprintf('DEEP3, LANE1 sample sheet save to %s', outpath))
write_delim(deep3_low2_lane3_sampleSheet, paste0(outpath, 'deep3_low2_lane3_sampleSheet.txt'), delim = '\t')

# Amount of data required by city
output_gb <- 800

# Subtract Parents and F2 from number of samples
num_samples <- nrow(deep3_low2_lane3_sampleSheet) - nrow(deep3_low2_lane3_sampleSheet %>% filter(is.na(city)))

deep3_low2_lane3_sampleSheet %>% 
  group_by(lib_split) %>% 
  summarise(count = n()) %>% 
  mutate(data_amount = case_when(lib_split == 'Parent_HomDel' ~ count * 10,
                                 lib_split == 'Parent_HomDom' ~ count * 10,
                                 lib_split == 'F2' ~ count * 5,
                                 TRUE ~ ((output_gb - 65) / num_samples) * count))  # Subtract coverage going to Parents and F2

