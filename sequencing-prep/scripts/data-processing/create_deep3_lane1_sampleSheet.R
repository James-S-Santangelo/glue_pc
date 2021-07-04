# Script to create sample sheet for plants sequenced as part of DEEP3, LANE1
# Removing samples with low concentrations and those with no more library
# We had an issue with the first attempt at pooling these libraries (many failed QC)
# so we had to repeat the pooling. Also including the index sequences now, which will
# facilitate filling out the Sample Information Form for the sequencing centre. 

# Load in data and add native vs. introduced range
isdup <- function (x) duplicated (x)
deep3_lane1_sampleSheet <- read_csv('data/clean/deep3/deep3_lane1_dilutions.csv') %>% 
  
  dplyr::select(continent, city, pop, individual, site, plantID, lane, i5_primer_index_seq, i7_primer_index_seq) %>% 
  filter(lane == 1) %>% 
  rename('sample' = 'plantID') %>% 
  
  # Add "_dup" suffix to duplicated samples (second in pair only).
  # This was done manually when creating the SIF. Only 2 duplicate samples
  mutate(is_dup = isdup(sample),
         sample = ifelse(is_dup, paste0(sample, '_dup'), sample)) %>% 
  dplyr::select(-is_dup) %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                           city == 'Tehran' ~ 'Native',
                           TRUE ~ 'Introduced'))
  
outpath <- 'resources/'
print(sprintf('DEEP3, LANE1 sample sheet save to %s', outpath))
write_delim(deep3_lane1_sampleSheet, paste0(outpath, 'deep3_lane1_sampleSheet.txt'), delim = '\t')

