# Script used to extract Human Influence Index

# Load all population-mean dataframes as list
inpath <- 'data/clean/popMeans_allCities/'
popMeans_dfList <- create_df_list(inpath)

#' Extracts population-specific Human Influence Index (HII) from raster
#' 
#' V2 of HII available from https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic
#' 
#' @param df Dataframe with population latitude and longitude coordinates
#' @param raster_path Character string containing full path to Raster file
#' 
#' @return dataframe with `hii` as column
extract_hii <- function(df, raster_path, outpath){
  
  # Get city name
  city <- df %>% pull(city) %>% unique()
  
  # Add fake coordinates for missing data
  df <- df %>% 
    mutate(population_longitude = ifelse(is.na(population_longitude), -999, population_longitude),
           population_latitude = ifelse(is.na(population_latitude), -999, population_latitude))
  
  # Load in raster for country
  raster <- raster::raster(raster_path)
  
  # Create spatial point dataframe from latitude and longitude
  spdf <- SpatialPointsDataFrame(coords = df %>% 
                                   dplyr::select(population_longitude, population_latitude), 
                                 proj4string = raster@crs, 
                                 data = df)
  
  # Extract GMIS data for population
  hii_data <- raster::extract(x = raster, y = spdf, method = 'simple')
  
  # Add column with HII values
  df_out <- df %>% 
    mutate(hii = hii_data,
           hii = ifelse(hii == 255, NA, hii))
  
  write_csv(df_out, paste0(outpath, city, '_hii.csv'))
}

# Extract HII for all cities and population
raster_path <- '~/Downloads/hfp-global-geo-grid/hfp_global_geo_grid/hf_v2geo/w001001.adf'
outpath <- 'data/raw/environmental_data/hii/'
dir.create(outpath)
purrr::walk(popMeans_dfList, extract_hii, raster_path = raster_path, outpath = outpath)
