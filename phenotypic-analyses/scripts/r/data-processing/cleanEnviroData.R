# Script to clean/filter all environmental variables collected for each city. 
# Required prior to merging environmental data with population-mean HCN frequency datasets
#
# Author: James. S. Santangelo

#### ANNUAL AI ####

inpath = "data/raw/environmental_data/Extracted_annual_Aridity_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/annualAI/"
purrr::walk(df_list, filter_environmental_data, variable = "annualAI", 
            outpath = outpath)

#### ANNUAL PET ####

inpath = "data/raw/environmental_data/Extracted_annual_PET_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/annualPET/"
purrr::walk(df_list, filter_environmental_data, variable = "annualPET", 
            outpath = outpath)

#### DEM ####

inpath = "data/raw/environmental_data/Extracted_DEM_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/DEM/"
purrr::walk(df_list, filter_environmental_data, variable = "DEM", 
            outpath = outpath)

#### GMIS ####

inpath = "data/raw/environmental_data/Extracted_GMIS_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/GMIS/"
dir.create(outpath)
purrr::walk(df_list, filter_environmental_data, variable = "GMIS", 
            outpath = outpath)

#### NDSI ####

inpath = "data/raw/environmental_data/Extracted_NDSI_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/NDSI/"
purrr::walk(df_list, filter_environmental_data, variable = "NDSI", 
            outpath = outpath)

#### SUMMER LST ####

inpath = "data/raw/environmental_data/Extracted_summer_LST_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/summerLST/"
purrr::walk(df_list, filter_environmental_data, variable = "summerLST", 
            outpath = outpath)

#### SUMMER NDVI ####

inpath = "data/raw/environmental_data/Extracted_summer_NDVI_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/summerNDVI/"
purrr::walk(df_list, filter_environmental_data, variable = "summerNDVI", 
            outpath = outpath)

#### WINTER LST ####

inpath = "data/raw/environmental_data/Extracted_winter_LST_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/winterLST/"
purrr::walk(df_list, filter_environmental_data, variable = "winterLST", 
            outpath = outpath)

#### WINTER NDVI ####

inpath = "data/raw/environmental_data/Extracted_winter_NDVI_filtered/"
df_list <- create_df_list(inpath)
outpath = "data/clean/environmental_data/winterNDVI/"
purrr::walk(df_list, filter_environmental_data, variable = "winterNDVI", 
            outpath = outpath)

