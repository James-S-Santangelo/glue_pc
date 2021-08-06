# Script to clean post-PCR library concentrations for the third lane of DEEP3 sequencing
#
# Author: James S. Santangelo

# Read in data with library concentrations and other data
# Note this differ slightly from the one initially created since
# some LOW2 cities samples were manually removed to make room for 
# an additional 10 plants of known cyanotype. 
# These 10 plants will be sequenced here as part of Lane 3
prepped_library_df <- read_csv('data/raw/20210806_plantsToPrep_deep3_low2_thirdLanePrepped.csv') %>% 
  
  # Select relevant columns
  dplyr::select(group, continent, city, pop, individual, site, plantID, max_qubit, leftover, Bioruptor_label,
                lane, i5, i7, 'Date prepped', Notes, Qubit_toUse) %>% 
  
  # Rename columns
  rename('date_prepped' = 'Date prepped') %>% 
  
  # Include only lane 2
  filter(lane == 3)

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
  mutate(plantID = str_replace(plantID, '^Punta_Arenas', 'Punt_Ar'),
         plantID = str_replace(plantID, '^Thessaloniki', 'Thessa')) %>% 
  
  # Create column used for splitting libraries.
  # Single library for each of Marc's parent plants (N = 2)
  # One libraries with all 8 F2s (N = 1)
  # One library for each city in deep3 (N = 8)
  # One library with the 19 cities from low2 (N = 1)
  # Total of 12 libraries
  mutate(lib_split = case_when(plantID == 'Parent_HomDom' ~ 'Parent_HomDom',
                               plantID == 'Parent_HomDel' ~ 'Parent_HomDel',
                               str_detect(plantID, '^F2_') ~ 'F2',
                               group == 'deep3' ~ city,
                               group == 'low2' ~ 'low2',
                               TRUE ~ NA_character_)) %>% 
  filter(!(is.na(lib_split)))

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
write_csv(prepped_library_df_mod, paste0(outpath, 'deep3_low2_lane3_libraryConcentrations.csv'),
          col_names = TRUE)  

