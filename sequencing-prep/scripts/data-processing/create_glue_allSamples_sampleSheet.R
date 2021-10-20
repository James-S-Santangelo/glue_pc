# Script to combine sample sheets from all sequencing runs
#
# Autor: James S. Santangelo

# Load in sample sheets
low1_ss <- read_delim('resources/low1_sampleSheet.txt', delim = '\t') %>% 
  filter(!(city == 'Toronto')) %>% 
  mutate(library = 'glue_low1',
         lane = 1)
cols = names(low1_ss)  # Keep all columns from LOW1
deep3_lane1_ss <- read_delim('resources/deep3_lane1_sampleSheet.txt', delim = '\t') %>% 
  mutate(library = 'deep3') %>% 
  select(cols)
deep3_lane2_ss <- read_delim('resources/deep3_lane2_sampleSheet.txt', delim = '\t') %>% 
  mutate(library = 'deep3') %>% 
  select(cols)
deep3_low2_lane3_ss <- read_delim('resources/deep3_low2_lane3_sampleSheet.txt', delim = '\t') %>% 
  filter(!is.na(individual)) %>%  # Removes samples not part of GLUE
  rename('library' = 'lib_split') %>% 
  mutate(library = ifelse(library == 'low2', 'low2', 'deep3')) %>% 
  select(cols)
toronto_gwsd <- read_delim('resources/sequencedPlants_phenotypesHabitat.txt', delim = '\t') %>% 
  rename('sample' = 'Sample',
         'site' = 'Habitat',
         'pop' = 'Population',
         'individual' = 'Plant') %>% 
  dplyr::select(-HCN_Result, -Locus.Li, -Locus.Ac) %>% 
  mutate(continent = 'NAM',
         city = 'Toronto',
         range = 'Introduced',
         library = 'Toronto',
         lane = 1)

# Combine sample sheets
sample_sheet <- bind_rows(low1_ss, deep3_lane1_ss, deep3_lane2_ss, deep3_low2_lane3_ss, toronto_gwsd) %>% 
  filter(!(sample == 'Medellin_39_10_dup'))  ## This sample was missing from sequencing

write_delim(sample_sheet, 'resources/glue_allSamples_sampleSheet.txt', delim = '\t')


