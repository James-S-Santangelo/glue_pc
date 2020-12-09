# Script to merge indices with samples

# Load libraries
library(tidyverse)

# Load data with samples
lane1 <- read_csv("resources/illumina-sequencing/library-preps/pooling/GWSD_Toronto_Lane1.csv")
lane2 <- read_csv("resources/illumina-sequencing/library-preps/pooling/GWSD_Toronto_Lane2.csv")

# Load dataframes with indices
i5_indices <- read_csv("resources/illumina-sequencing/reference/indices/iTru5_forward-indices.csv") %>% 
  rename("i5_index" = "Primer_Name") %>% 
  select(i5_index, Primer_Index_Sequence)
i7_indices <- read_csv("resources/illumina-sequencing/reference/indices/iTru7_reverse-indices.csv") %>% 
  rename("i7_index" = "Primer_Name") %>% 
  select(i7_index, Primer_Index_Sequence)

#### FUNCTIONS ####
merge_indices <- function(df, i5_indices, i7_indices){
  
  
  df_out <- df %>% 
    select(City, Population, Plant, plant_id,
           i5_index, i7_index, Column, Row, library_volume_uL, TE_vol_uL) %>% 
    left_join(i5_indices, by = "i5_index") %>%
    rename("i5_seq" = "Primer_Index_Sequence") %>% 
    left_join(i7_indices, by = "i7_index") %>% 
    rename("i7_seq" = "Primer_Index_Sequence")
  
  return(df_out)
}

lane1_out <- merge_indices(lane1, i5_indices, i7_indices)
lane2_out <- merge_indices(lane2, i5_indices, i7_indices)

write_csv(lane1_out, path = "resources/illumina-sequencing/library-preps/pooling/GWSD_Toronto_lane1_withIndices.csv")
write_csv(lane2_out, path = "resources/illumina-sequencing/library-preps/pooling/GWSD_Toronto_lane2_withIndices.csv")

