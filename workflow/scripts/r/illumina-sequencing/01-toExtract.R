# Script to determine from which Toronto populations DNA should be extracted
#
# Author: James S. Santangelo
# Date: July 26, 2019

# Load in required packages
library(tidyverse)

#Load in dataset containing Toronto plants
Tor_plants <- read_csv("resources/illumina-sequencing/reference/allPlants_Toronto.csv") %>% 
  mutate(Transect = case_when(Transect == "A" ~ "North",
                              Transect == "B" ~ "West",
                              Transect == "C" ~ "East"))

# Get 7 most urban populations
Urban_pops <- Tor_plants %>% 
  arrange(Distance, .by_group = TRUE) %>%
  distinct(Distance, Population, Transect, Lat.pop, Long.pop) %>%
  top_n(-7, Distance) %>% 
  mutate(Habitat = "Urban")

# Get 7 most rural populations by transect
Rural_pops <- Tor_plants %>% 
  group_by(Transect) %>% 
  arrange(desc(Distance), .by_group = TRUE) %>%
  distinct(Distance, Population, Transect, Lat.pop, Long.pop) %>%
  top_n(7, Distance) %>% 
  mutate(Habitat = "Rural")
  
# Get seven suburban populations
# Pops are 3 km on either side of transect midpoint
Suburban_pops <- Tor_plants %>% 
  group_by(Transect) %>% 
  mutate(Mid = (max(Distance) - min(Distance)) / 2,
         is_suburban = ifelse(Distance < Mid + 3 & Distance > Mid - 3, 1, 0)) %>% 
  filter(is_suburban == 1) %>% 
  distinct(Distance, Population, Transect, Lat.pop, Long.pop) %>% 
  # Transect B had 7 populations. Remove closest and furthest populations from city centre
  filter(!(Transect == "West" & (Distance == min(Distance) | Distance == max(Distance)))) %>% 
  mutate(Habitat = "Suburban")

all_tor_pops <- bind_rows(Urban_pops, Rural_pops, Suburban_pops) %>% 
  rename("latitude" = "Lat.pop",
         "longitude" = "Long.pop")

# Write toronto populations to disk
write_csv(all_tor_pops, "resources/illumina-sequencing/library-preps//01_torontoPops_toExtract.csv")
