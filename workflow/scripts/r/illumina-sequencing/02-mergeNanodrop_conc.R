# Script merge Nanodrop concentrations and purity with extraction data (+ Qubit)
#
# Author: James S. Santangelo
# Date: Sept. 19, 2019

# Load required packages
library(tidyverse)

# Load in extraction data with qubit readings
extraction_data <- read_csv("data/illumina-sequencing/raw/allPlants_Toronto_extractionData.csv") %>% 
  select(City:Qubit_conc)


#' Creates list with dataframes as elements for all CSVs in inpath
#' @param inpath Path to directory containing CSVs to load
#' @return df_list. List of dataframes, each from a different CSV
create_df_list <- function(inpath){
  
  # Get all csv files in inpath
  files <- dir(inpath, pattern = "*.csv")
  
  # read in all the files, appending the path before the filename
  df_list <- files %>%
    map(~ read_csv(file.path(inpath, .), skip = 2))
  
  return(df_list)
}

clean_nanodrop_data <- function(df){
  
  df_mod <- df %>% 
    select(Name, "260/280", "ng/µL") %>% 
    rename("260_280" = "260/280", "nanodrop_conc" = "ng/µL") %>% 
    filter(!is.na(nanodrop_conc)) %>% 
    mutate(Name = str_replace_all(Name, "Tor", "")) %>% 
    separate(., Name, into = c("Population", "Plant"), sep = "[.]") %>% 
    mutate(Population = as.numeric(Population),
           Plant = as.numeric(Plant))
  
  return(df_mod)
}

# Read in all Nanodrop dfs
nanodrop_dfs <- create_df_list("data/illumina-sequencing/reference/nanodrop/")

# Clean all Nanodrop dfs
nanodrop_dfs_mod <- purrr::map(nanodrop_dfs, clean_nanodrop_data)

# Merge nanoodrop dfs into single df
nanodrop_merged <- bind_rows(nanodrop_dfs_mod)

extraction_data_merged <- extraction_data %>% 
  left_join(nanodrop_merged, by = c("Population", "Plant"))

write_csv(extraction_data_merged, "data/illumina-sequencing/library-preps/02_allPlants_Toronto_extractionData_allQC.csv")
