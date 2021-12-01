# Script to add population-level environmental data to population-mean HCN frequency datasets
#
# Author: James S. Santangelo

# Get list of population-mean HCN frequency datasets.
inpath <- "data/clean/popMeans_allCities/"
popMeans_dfList <- create_df_list(inpath)

# Write population mean datasets with enviro data to disk
# Additionally writes text file with missing data by city (if any)
outpath = "data/clean/popMeans_allCities_withEnviro/"
err_file = paste0(outpath, "missingEnviroData.txt")
file.create(err_file)
purrr::walk(popMeans_dfList, add_enviro_data, 
            outpath = outpath, err_file = err_file)