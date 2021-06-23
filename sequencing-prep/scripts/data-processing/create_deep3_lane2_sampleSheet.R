# Script to create sample sheet for plants sequenced as part of DEEP3, LANE2

# Load in data and add native vs. introduced range
deep3_lane2_sampleSheet <- read_csv('data/clean/deep3/deep3_lane2_dilutions.csv') %>% 

  dplyr::select(continent, city, pop, individual, site, plantID, lane, i5_primer_index_seq, i7_primer_index_seq) %>% 
  rename('sample' = 'plantID') %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                           city == 'Tehran' ~ 'Native',
                           TRUE ~ 'Introduced'))
  

outpath <- 'resources/'
print(sprintf('DEEP3, LANE1 sample sheet save to %s', outpath))
write_delim(deep3_lane2_sampleSheet, paste0(outpath, 'deep3_lane2_sampleSheet.txt'), delim = '\t')
