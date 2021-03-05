# Script to create sample sheet for plants sequenced as part of LOW1

# 20 samples from Toronto
tor_samples <- read_csv('resources/toronto_samples_forGlue.csv')

# Load in data and add native vs. introduced range
low1_sampleSheet <- read_csv('data/clean/low1/low1_libraryConcentrations.csv') %>% 
  dplyr::select(continent, city, pop, individual, site, plantID) %>% 
  rename('sample' = 'plantID') %>% 
  mutate(range = case_when(continent == 'EU' ~ 'Native',
                          city == 'Tehran' ~ 'Native',
                          TRUE ~ 'Introduced')) %>% 
  bind_rows(., tor_samples)

outpath <- 'resources/'
dir.create(outpath, showWarnings = FALSE)
print(sprintf('LOW1 sample sheet save to %s', outpath))
write_delim(low1_df, paste0(outpath, 'low1_sampleSheet.txt'), delim = '\t')