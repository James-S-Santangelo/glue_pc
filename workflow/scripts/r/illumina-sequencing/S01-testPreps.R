# Script to select plants for test and optimization preps

# Load in data with extracted Toronto plants
torPlants <- read_csv("resources/illumina-sequencing/library-preps/03_allPlants_Toronto_allData.csv") %>% 
  filter(!is.na(Qubit_conc)) %>% 
  mutate(plant_id = paste(Population, Plant, sep = "_"))

# Load plants to library prep
allPlant_toPrep <- read_csv("resources/illumina-sequencing/library-preps/04_allPlants_toPrep.csv")

final_vol <- 50 # Volume required for shearing
final_conc <- 10 # Concentration required for shearing
initial_vol <- 40 # Assume I have 40 uL for all sampled (I should have a bit more)

# Need at least this concentration for samples
initial_conc <- (final_vol * final_conc) / initial_vol

# 8 Plant on which to test multichannel library prep
# These plant have been completely exhausted while optimizing the shearing
# The library prep was ruined due to insufficient enzymes 
set.seed(44)
testPrep <- torPlants %>% 
  filter(Population %in% c(6, 7, 42, 43) &
           Qubit_conc > initial_conc &
           !(plant_id %in% allPlant_toPrep$plant_id)) %>% 
  group_by(Population) %>% 
  sample_n(2) %>% 
  arrange(Population, Plant) %>% 
  select(City, Population, Plant, plant_id, Plate, Qubit_conc, '260_280') %>% 
  mutate(initial_vol = round((5 * final_conc) / Qubit_conc, 2),
         TE_vol = round(5 - initial_vol, 2))

write_csv(testPrep, path = "resources/illumina-sequencing/library-preps/optimizing-testing/01_testPrep_plants.csv", col_names = TRUE)

#### TEST SHEARING ####

# Select plants used for optimizing the shearing process
# Shearing test successful but plants were not prepped

# Set final concentrations for high concentration and low concentration treatments
final_vol_BR_course <- 25
final_conc_BR_course_low <- 10
final_conc_BR_course_high <- 20
intial_vol_BR_course <- 40

# Required initial concentration of samples
initial_conc_BR_course_low <- (final_vol_BR_course*final_conc_BR_course_low) / (intial_vol_BR_course)
initial_conc_BR_course_high <- (final_vol_BR_course*final_conc_BR_course_high) / (intial_vol_BR_course)

# Set seed
set.seed(45)

# Select plants
bioruptor_timeCourse <- torPlants %>% 
  filter(Qubit_conc > initial_conc_BR_course_high &
           !(plant_id %in% allPlant_toPrep$plant_id) & # Plant not being used in Pilot
           !(plant_id %in% testPrep$plant_id) & # Plant not already test sheared
           `260_280` > 1.8 & `260_280` < 2 & # Only good quality DNA
           Plate == 10 & # All from single plate with many good plants
           Qubit_conc <= 105) %>% # Get rid of really high concentration plants
  select(City, Population, Plant, plant_id, Plate, Qubit_conc, `260_280`) %>%
  sample_frac(1) %>% # Randomly reshuffle plants
  
  # Create treatment for concentration based on initial concentrations of available samples
  mutate(shear_conc = case_when(Qubit_conc > 20 ~ final_conc_BR_course_high,
                                Qubit_conc <= 20 ~ final_conc_BR_course_low),
         cycles = c(4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 4, 5), # Number of shearing cycles
         on = c(5, 10, 15, 5, 10, 15, 5, 10, 15, 5, 10, 15), # Time ON per cycle
         
         # Get required initial sample volume based on sample's treatment
         initial_vol = case_when(shear_conc == final_conc_BR_course_low ~ round((25 * final_conc_BR_course_low) / Qubit_conc, 2),
                                 shear_conc == final_conc_BR_course_high ~ round((25 * final_conc_BR_course_high) / Qubit_conc, 2)),
         
         # Amount of TE to be added to final vol. of 25 
         TE_vol = round(final_vol_BR_course - initial_vol, 2)) %>% 
  mutate(treatment = paste(shear_conc, cycles, on, sep = '_')) %>% 
  as.data.frame() %>% 
  arrange(Population, Plant)

write_csv(bioruptor_timeCourse, path = 'resources/illumina-sequencing/library-preps/optimizing-testing/bioruptor_timeCourse.csv', col_names = TRUE)

#### TEST PREP: ROUND 2 ####

set.seed(44)
testPlants_02 <- torPlants %>% 
  filter(Qubit_conc > initial_conc &
           !(plant_id %in% allPlant_toPrep$plant_id) & # Plant not being used in Pilot
           !(plant_id %in% testPrep$plant_id) & # Plant not already test sheared
           !(plant_id %in% bioruptor_timeCourse$plant_id) &
           # `260_280` > 1.8 & `260_280` < 2 & # Only good quality DNA
           Plate == 10 & # All from single plate with many good plants
           Qubit_conc <= 150) %>% # Get rid of really high concentration plants
  select(City, Population, Plant, plant_id, Plate, Qubit_conc, `260_280`) %>%
  sample_n(8) %>% 
  mutate(cycles = c(4, 4, 4, 4, 5, 5, 5, 5)) %>% 
  arrange(Population, Plant) %>% 
  mutate(initial_vol = round((final_vol * final_conc) / Qubit_conc, 1),
         TE_vol = round(final_vol - initial_vol, 1),
         tube = 1:8) 

write_csv(testPlants_02, path = 'resources/illumina-sequencing/library-preps/optimizing-testing/02-03_testPrep_plants.csv', col_names = TRUE)

#### TEST PREP: ROUND 3 ####

# This prep used the same plants from 02_test prep and are thus not
# selected here

#### TEST PREP: ROUND 4 ####

allPlant_toPrep %>% 
  ungroup() %>% 
  group_by(Plate) %>% 
  tally()

# Will use plants from plate 3 since it has the fewest plants being
# unsed in the actual pilot (with the exception of plate 10 which has
# since been mostly used for trials)

# Will be testing the efficacy of 3 vs. 4 bioruptor cycles plus
# whether size selection could be done after fragmentation or should
# be left to after PCR

testPlants_04 <- torPlants %>% 
  filter(Qubit_conc > initial_conc &
           !(plant_id %in% allPlant_toPrep$plant_id) & # Plant not being used in Pilot
           !(plant_id %in% testPrep$plant_id) & # Plant not already test sheared
           !(plant_id %in% bioruptor_timeCourse$plant_id) & # Plant not used in initial bioruptor test
           !(plant_id %in% testPlants_02$plant_id)) %>% # Plant not used in second bioruptor test
  filter(Plate == 3) %>% # Will use plant from plate 3
  select(City, Population, Plant, plant_id, Plate, Qubit_conc, `260_280`) %>% 
  top_n(2, Qubit_conc) %>% 
  mutate(initial_vol = round((25 * final_conc) / Qubit_conc, 1),
         TE_vol = round(25 - initial_vol, 1)) 

write_csv(testPlants_04, path = 'resources/illumina-sequencing/library-preps/optimizing-testing/04_testPrep_plants.csv', col_names = TRUE)


#### TEST PREP: ROUND 5 ####

allPlant_toPrep %>% 
  ungroup() %>% 
  group_by(Plate) %>% 
  tally()

# Will use plants from plate 3 since it has the fewest plants being
# unsed in the actual pilot (with the exception of plate 10 which has
# since been mostly used for trials)

# These plants will be used by Sophie to test the library prep protocol.

SKtestPlants <- torPlants %>% 
  filter(Qubit_conc > initial_conc &
           !(plant_id %in% allPlant_toPrep$plant_id) & # Plant not being used in Pilot
           !(plant_id %in% testPrep$plant_id) & # Plant not already test sheared
           !(plant_id %in% bioruptor_timeCourse$plant_id) & # Plant not used in initial bioruptor test
           !(plant_id %in% testPlants_02$plant_id) &
           !(plant_id %in% testPlants_04$plant_id)) %>% # Plant not used in second bioruptor test
  filter(Plate == 3) %>% # Will use plant from plate 3
  select(City, Population, Plant, plant_id, Plate, Qubit_conc, `260_280`) %>% 
  top_n(4, Qubit_conc) %>%
  mutate(initial_vol = round((30 * final_conc) / Qubit_conc, 1),
         TE_vol = round(30 - initial_vol, 1)) 

write_csv(SKtestPlants, path = 'resources/illumina-sequencing/library-preps/optimizing-testing/05_SK_testPrep_plants.csv', 
          col_names = TRUE)
