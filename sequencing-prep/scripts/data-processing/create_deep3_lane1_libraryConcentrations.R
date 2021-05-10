# Script to clean post-PCR library concentrations for the first lane of DEEP3 sequencing
#
# Author: James S. Santangelo

# Read in data with library concentrations and other data
# Note there are two sheets for this in the raw data. Only the one listed here should be used. 
# Only the Qubit3_HS_reEst column is needed. Previous Qubit estimates were incorrect and had to
# be re-done.
prepped_library_df <- read_csv('data/raw/20210510_plantsToPrep_deep3_low2_firstLanePrepped_qubitReEst.csv') %>% 
  
  # Select relevant columns
  dplyr::select(continent, city, pop, individual, site, plantID, max_qubit, leftover, Bioruptor_label,
         lane, i5, i7, 'Date prepped', Notes, contains('HS')) %>% 
  
  # Rename columns
  rename('date_prepped' = 'Date prepped',
         'qubit1_hs' = 'Qubit1_HS (ng/ul)',
         'qubit2_hs' = 'Qubit2_HS (ng/ul)',
         'qubit3_hs_reEst' = 'Qubit3_HS_reEst') %>% 
  
  # Include only lane 1
  filter(lane == 1) %>% 
  
  # Get max qubit for prepped libraries
  # In this case, this is only done for Albuquerque since previous pooling was successful
  # For other cities, only the 'qubit3_hs_reEst' column is used. 
  mutate(library_qubit = case_when(city == 'Albuquerque' ~ pmax(qubit1_hs, qubit2_hs, na.rm = TRUE),
                                  TRUE ~ qubit3_hs_reEst))

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
  duplicated() %>% 
  sum()

outpath <- 'data/clean/deep3/'
print(sprintf('DEEP3 with library concentrations saved to %s', outpath))
write_csv(prepped_library_df_mod, paste0(outpath, 'deep3_lane1_libraryConcentrations.csv'),
          col_names = TRUE)  


test <- prepped_library_df %>% filter(!(city == 'Albuquerque')) %>% 
  mutate(previous_qubit = pmax(qubit1_hs, qubit2_hs, na.rm = TRUE))
cor(test$previous_qubit, test$library_qubit)
plot(test$previous_qubit, test$library_qubit)
