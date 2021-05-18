# Script to create sample sheet for plants sequenced as part of DEEP3, LANE1
# Removing samples with low concentrations and those with no more library
# We had an issue with the first attempt at pooling these libraries (many failed QC)
# so we had to repeat the pooling. Also including the index sequences now, which will
# facilitate filling out the Sample Information Form for the sequencing centre. 

# Load in data and add native vs. introduced range
deep3_lane1_sampleSheet <- read_csv('data/clean/deep3/deep3_lane1_dilutions.csv') %>% 
  
  # Remove samples that were not included in the end
  # Low concentration samples from Antwerp
  filter(!(city == 'Antwerp' & is.na(library_volume_uL))) %>% 
  # Low concentration samples from Armidale
  filter(!(city == 'Armidale' & is.na(library_volume_uL))) %>% 
  # Remove 2 libraries with no remaining DNA
  filter(!(Bioruptor_label %in% c('ARM-71','ARM-72'))) %>% 
  dplyr::select(continent, city, pop, individual, site, plantID, lane, i5_primer_index_seq, i7_primer_index_seq) %>% 
  filter(lane == 1) %>% 
  rename('sample' = 'plantID') %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                           city == 'Tehran' ~ 'Native',
                           TRUE ~ 'Introduced'))
  

outpath <- 'resources/'
print(sprintf('DEEP3, LANE1 sample sheet save to %s', outpath))
write_delim(deep3_lane1_sampleSheet, paste0(outpath, 'deep3_lane1_sampleSheet.txt'), delim = '\t')
