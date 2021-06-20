# Script to clean post-PCR library concentrations for the second lane of DEEP3 sequencing
#
# Author: James S. Santangelo

# Read in data with library concentrations and other data
prepped_library_df <- read_csv('data/raw/20210620_plantsToPrep_deep3_low2_secondLanePrepped.csv') %>% 
  
  # Select relevant columns
  dplyr::select(continent, city, pop, individual, site, plantID, max_qubit, leftover, Bioruptor_label,
         lane, i5, i7, 'Date prepped', Notes, contains('HS')) %>% 
  
  # Rename columns
  rename('date_prepped' = 'Date prepped',
         'qubit1_hs' = 'Qubit1_HS (ng/ul)',
         'qubit2_hs' = 'Qubit2_HS (ng/ul)',
         'qubit3_hs_reEst' = 'Qubit3_HS_reEst') %>% 
  
  # Include only lane 1
  filter(lane == 2) %>% 
  
  # Get max qubit for prepped libraries
  mutate(library_qubit = pmax(qubit1_hs, qubit2_hs, qubit3_hs_reEst, na.rm = TRUE))

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
isdup <- function (x) duplicated (x) | duplicated (x, fromLast = TRUE)
prepped_library_df_mod %>% 
  mutate(primer_insex_seq_comb = paste0(i5_primer_index_seq, '_', i7_primer_index_seq)) %>% 
  dplyr::select(primer_insex_seq_comb, Bioruptor_label, plantID) %>% 
  mutate(is_dup = isdup(primer_insex_seq_comb)) %>% 
  filter(is_dup) %>% 
  View()

outpath <- 'data/clean/deep3/'
print(sprintf('DEEP3 with library concentrations saved to %s', outpath))
write_csv(prepped_library_df_mod, paste0(outpath, 'deep3_lane2_libraryConcentrations.csv'),
          col_names = TRUE)  

