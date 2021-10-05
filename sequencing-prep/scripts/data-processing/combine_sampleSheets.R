# Script to combine sample sheets from all sequencing runs
#
# Autor: James S. Santangelo

# Load in sample sheets
low1_ss <- read_delim('resources/low1_sampleSheet.txt', delim = '\t') %>% 
  mutate(library = ifelse(city == 'Toronto', 'Toronto', 'glue_low1'),
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

# Combine sample sheets
sample_sheet <- bind_rows(low1_ss, deep3_lane1_ss, deep3_lane2_ss, deep3_low2_lane3_ss)



