# Script to generate population-mean HCN datasets for all GLUE cities
#
# Author: James S. Santangelo

# Create list with paths to dataframes
inpath <- "data/clean/individualPlant_allCities/"
df_list <- create_df_list(inpath)

# Generate population mean datasets and store all as list
city_centers <- read_csv("data/clean/latLong_cityCenters_clean.csv")
popMeans_dfList <- purrr::map(df_list, generate_popMeans, 
                              city_centers = city_centers)

# Write pop means without environmental data. Used for pulling
# environmental data using Alex Tong's custom Python scripts
outpath = "data/clean/popMeans_allCities/"
purrr::walk(df_list, generate_popMeans, 
           city_centers = city_centers,
           outpath = outpath)


#### ADD MARC'S POP MEANS DATA ####

# recreate list with paths to MTJJ dataframes
inpath <- "data/raw/mtjj_jss/mtjj_popMeans/"
df_list <- create_df_list(inpath)

plantsPerPop <- read_csv('data/raw/johnson_2018_plantsPerPop.csv')
mtjj_popMeans <- purrr::map(df_list, clean_mtjj_popMeans, 
                             city_centers = city_centers,
                            plantsPerPop = plantsPerPop)

# Write pop means without environmental data. Used for pulling
# environmental data using Alex Tong's custom Python scripts
outpath = "data/clean/popMeans_allCities/"
purrr::walk(df_list, clean_mtjj_popMeans, 
            city_centers = city_centers,
            plantsPerPop = plantsPerPop,
            outpath = outpath)
