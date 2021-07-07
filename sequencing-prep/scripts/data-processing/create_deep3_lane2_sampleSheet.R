# Script to create sample sheet for plants sequenced as part of DEEP3, LANE2

# Load in data and add native vs. introduced range
isdup <- function (x) duplicated (x)
deep3_lane2_sampleSheet <- read_csv('data/clean/deep3/deep3_lane2_dilutions.csv') %>% 

  dplyr::select(continent, city, pop, individual, site, plantID, lane, i5_primer_index_seq, i7_primer_index_seq) %>% 
  rename('sample' = 'plantID') %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                           city == 'Tehran' ~ 'Native',
                           TRUE ~ 'Introduced')) %>% 
  
  # Shorten Palmerston North sample IDs since they're too long for Novogene Sample Sheet
  mutate(sample = str_replace(sample, '^Palmerston_North', 'Palm_North')) %>% 
  
  # Add '_dup_' suffix to duplicated samples
  mutate(is_dup = isdup(sample),
         sample = ifelse(is_dup, paste0(sample, '_dup'), sample)) %>% 
  dplyr::select(-is_dup)
  

outpath <- 'resources/'
print(sprintf('DEEP3, LANE1 sample sheet save to %s', outpath))
write_delim(deep3_lane2_sampleSheet, paste0(outpath, 'deep3_lane2_sampleSheet.txt'), delim = '\t')

# Amount of data required by city
output_gb <- 800
num_samples <- nrow(deep3_lane2_sampleSheet)
deep3_lane2_sampleSheet %>% 
  group_by(city) %>% 
  summarise(count = n()) %>% 
  mutate(data_amount = (output_gb / num_samples) * count)
