# Script to concatenate DEEP-3 and LOW-2 and assign lanes

# Load datasets
deep3 <- read_csv('data/clean/deep3/plantsToPrep_deep3.csv')
low2 <- read_csv('data/clean/low2/plantsToPrep_low2.csv')

# Concatenate dataframe
deep3_low2 <- bind_rows('deep3' = deep3, 'low2' = low2, .id = 'group')

# Figure out how many plants per lane
num_plants <- nrow(deep3_low2)
num_lanes <- 3
plants_per_lane <- num_plants / num_lanes

# Add lane breakdown
deep3_low2 <- deep3_low2 %>% 
  mutate(row_id = row_number(),
         lane = case_when(row_id <= plants_per_lane ~ 1,
                          row_id > plants_per_lane & row_id <= plants_per_lane * 2 ~ 2,
                          row_id > plants_per_lane * 2 ~ 3)) %>% 
  select(-row_id)

# Confirm number of plants in each lane
deep3_low2 %>% 
  group_by(lane) %>% 
  tally()

outpath <- 'data/clean/'
dir.create(outpath, showWarnings = FALSE)
print(sprintf('DEEP3 and LOW2 concatenated file with lanes save to %s', outpath))
write_csv(deep3_low2, paste0(outpath, 'plantsToPrep_deep3_low2_withLanes.csv'))
