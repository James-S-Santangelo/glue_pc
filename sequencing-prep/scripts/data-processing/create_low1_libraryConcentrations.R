# Script to clean CSV with prepped libraries prior to calculating dilution volumes for pooling

# Read in data with library concentrations and other data
prepped_library_df <- read_csv('data/raw/20210105_plantsToPrep_Low1_Low2_firstLanePrepped.csv') %>% 
  
  # Select relevant columns
  dplyr::select(continent, city, pop, individual, site, plantID, max_qubit, leftover, Bioruptor_label,
         'Batch/lane', i5, i7, 'Date prepped', Notes, ends_with('(HS)')) %>% 
  
  # Rename columns
  rename('lane' = 'Batch/lane',
         'date_prepped' = 'Date prepped',
         'qubit1_hs' = 'Qubit1_ng/ul (HS)',
         'qubit2_hs' = 'Qubit2_ng/ul(HS)',
         'qubit3_hs' = 'Qubit3_ng/ul(HS)',
         'qubit4_hs' = 'Qubit4_ng/ul(HS)') %>% 
  
  # Include only lane 1
  filter(lane == 1) %>% 
  
  # Get max qubit for prepped libraries
  mutate(max_library_qubit = pmax(qubit1_hs, qubit2_hs, qubit3_hs, qubit4_hs, na.rm = TRUE))

# Load in index sequences
i5_indices <- read_csv('resources/illumina_sequencing/iTru5_forward-indices.csv') %>% 
  dplyr::select(-X4, -X5, -X6, -Primer_Sequence)  # Remove extra columns
i7_indices <- read_csv('resources/illumina_sequencing/iTru7_reverse-indices.csv') %>% 
  dplyr::select(-Primer_Sequence)

# Join i5 and i7 primer idex sequences to dataframe
prepped_library_df_mod <- prepped_library_df %>% 
  left_join(., i5_indices, by = c('i5' = 'Primer_Name')) %>% 
  rename('i5_primer_index_seq' = 'Primer_Index_Sequence') %>%   
  left_join(., i7_indices, by = c('i7' = 'Primer_Name')) %>% 
  rename('i7_primer_index_seq' = 'Primer_Index_Sequence') %>% 
  mutate(plantID = str_replace(plantID, '^Buenos_Aires', 'Buen_Air'),
         plantID = str_replace(plantID, '^Christchurch', 'Chrchurch'),
         plantID = str_replace(plantID, '^Thessaloniki', 'Thessa'))

# Confirm than none of the index combinations are duplicated
prepped_library_df_mod %>% 
  mutate(primer_insex_seq_comb = paste0(i5_primer_index_seq, '_', i7_primer_index_seq)) %>% 
  dplyr::select(primer_insex_seq_comb) %>% 
  duplicated()

outpath <- 'data/clean/low1/'
print(sprintf('LOW1 with library concentrations saved to %s', outpath))
write_csv(prepped_library_df_mod, paste0(outpath, 'low1_libraryConcentrations.csv'),
          col_names = TRUE)  
