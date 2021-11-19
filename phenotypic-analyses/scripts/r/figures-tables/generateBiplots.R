# Load required packages
library(tidyverse)
source("scripts/r/utilityFunctions.R")

# Theme used for plotting

ng1 = theme(
  aspect.ratio = 0.7,
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(color = "black", size = 1),
  axis.line.y = element_line(color = "black", size = 1),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(color = "black", size = 15),
  axis.title = element_text(color = "black", size = 1),
  axis.title.y = element_text(vjust = 2, face = "bold", size = 18),
  axis.title.x = element_text(vjust = 0.1, face = "bold", size = 18),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  legend.position = "none",
  legend.direction = "vertical",
  legend.text = element_text(size = 13),
  legend.key = element_rect(fill = "white"),
  legend.title = element_text(size = 15, face = "bold"),
  legend.key.size = unit(1.0, "cm")
)

#### UNIVARIATE CLINES ####

# reate list with paths to dataframes
inpath <- "data/clean/popMeans_allCities_withEnviro/"
df_list <- create_df_list(inpath)

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/cline_biplots/"
purrr::walk(df_list, create_Biplot,
            outpath = outpath)
df_list[1]
#### ENVIRONMENTAL VARIABLES ####

## ANNUAL AI ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/annualAI/"
purrr::walk(df_list, create_Biplot, response_var = "annualAI_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## ANNUAL PET ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/annualPET/"
purrr::walk(df_list, create_Biplot, response_var = "annualPET_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## DEM ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/DEM/"
purrr::walk(df_list, create_Biplot, response_var = "DEM_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## GMIS ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/GMIS/"
purrr::walk(df_list, create_Biplot, response_var = "GMIS_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## NDSI ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/NDSI/"
purrr::walk(df_list, create_Biplot, response_var = "NDSI_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## SUMMER LST ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/summerLST/"
purrr::walk(df_list, create_Biplot, response_var = "summerLST_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## SUMMER NDVI ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/summerNDVI/"
purrr::walk(df_list, create_Biplot, response_var = "summerNDVI_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## WINTER NDVI ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/winterNDVI/"
purrr::walk(df_list, create_Biplot, response_var = "winterNDVI_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)

## WINTER LST ##

# Apply write_popMeans to each dataframe in datafame list
outpath <- "analysis/figures/environmental_biplots/winterLST/"
purrr::walk(df_list, create_Biplot, response_var = "winterLST_Mean", 
            predictor_var = "std_distance", 
            outpath = outpath)
