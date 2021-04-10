# Script to extract Impervious Surface (GMIS) for all sampled populations
# Requires locally available GMIS datasets, which can be downloaded from
# https://sedac.ciesin.columbia.edu/data/set/ulandsat-gmis-v1/data-download
# and unzipped into a single directory. Accessed April 9, 2021.
#
# Author: James S. Santangelo

#' Retrieves three-letter country code from country name
#'     
#' @param df Population-mean dataframe containing "Country" column
#'     
#' @return Country code as character vector
get_country_code <- function(df){
  
  # Get country
  country <- df %>% pull(Country) %>% unique()

  # Convert country to country code
  country_code <- case_when(country == 'Argentina' ~ 'ARG',
                            country == 'Australia' ~ 'AUS',
                            country == 'Belgium' ~ 'BEL',
                            country == 'Brazil' ~ 'BRA',
                            country == 'Canada' ~ 'CAN',
                            country == 'Chile' ~ 'CHL',
                            country == 'China' ~ 'CHN',
                            country == 'Colombia' ~ 'COL',
                            country == 'Ecuador' ~ 'ECU',
                            country == 'Finland' ~ 'FIN',
                            country == 'France' ~ 'FRA',
                            country == 'Germany' ~ 'DEU',
                            country == 'Greece' ~ 'GRC',
                            country == 'Iran' ~ 'IRN',
                            country == 'Japan' ~ 'JPN',
                            country == 'Mexico' ~ 'MEX',
                            country == 'Netherlands' ~ 'NLD',
                            country == 'Norway' ~ 'NOR',
                            country == 'New Zealand' ~ 'NZL',
                            country == 'Poland' ~ 'POL',
                            country == 'Portugal' ~ 'PRT',
                            country == 'Sweden' ~ 'SWE',
                            country == 'Switzerland' ~ 'CHE',
                            country == 'UK' ~ 'GBR',
                            country == 'USA' ~ 'USA',
                            country == 'South Africa' ~ 'ZAF',
                            TRUE ~ 'NA')

  return(country_code)
  
}

#' Matches country code to raster dataset and loads raster
#'     
#' @param country_code Three-letter country code for which to load raster
#' @param city City for which data will be extracted
#' @param inpath Path containing GMIS raster datasets for all countries
#'     
#' @return RasterLayer object. See ?raster
load_raster <- function(country_code, city, inpath){
  
  # Regex pattern to match country's raster
  # Handles a few edge cases for countries with multiple rasters
  if(country_code == 'USA' & !(city %in% c('Fairbanks', 'Anchorage'))){
    pattern <- '^USAw1_.*.tif$'
  }else if(country_code == 'USA' & city %in% c('Fairbanks', 'Anchorage')){
    pattern <- '^USAW3_.*.tif$'  
  }else if(country_code == 'NZL'){
    pattern <- '^NZLe_.*.tif$'  
  }else{
    pattern <- sprintf('^%s_.*.tif$', country_code)      
  }

  # Path to raster
  raster_file <- list.files(inpath, pattern = pattern, recursive = TRUE)
  raster_path <- paste0(inpath, raster_file)
  # print(raster_path)
  
  # Load and return raster
  raster <- raster::raster(raster_path)
  return(raster)
  
}

#' Extracts GMIS data for population
#'     
#' @param df Population-mean dataframe containing "Country" column
#' @param inpath Path containing GMIS raster datasets for all countries
#' @param outpath Path to which dataframe with GMIS values should be written
#' @param buffer Radius of buffer around population to include in calculation of GMIS
#'     
#' @return Country code as character vector
extract_gmis <- function(df, inpath, outpath, buffer = 100){

  # Get city name
  city <- df %>% pull(city) %>% unique()
  print(city)
  
  # Do not execute if file exists
  outpath <- paste0(outpath, city, '_GMIS.csv')
  if(file.exists(outpath)){
    print(sprintf('GMIS data already exists for %s', city))
  # Otherwise extract data
  }else{
    
    df <- df %>% 
      mutate(population_longitude = ifelse(is.na(population_longitude), -999, population_longitude),
             population_latitude = ifelse(is.na(population_latitude), -999, population_latitude))
    
    # Retrieve country code
    country_code <- get_country_code(df)
    
    # Load in raster for country
    raster <- load_raster(country_code, city, inpath)
    
    # Create spatial point dataframe from latitude and longitude
    spdf <- SpatialPointsDataFrame(coords = df %>% 
                                     dplyr::select(population_longitude, population_latitude), 
                                   proj4string = raster@crs, 
                                   data = df)
    
    # Extract GMIS data for population
    gmis_data <- raster::extract(x = raster, y = spdf, method = 'simple', buffer = buffer)
    gmis_data <- ifelse(is.na(gmis_data), 255, gmis_data)
    
    # Take mean GMIS across all cells included within buffer
    gmis_data_mod <- gmis_data %>% 
      # Convert GMIS values of 200 to 0, as per documentation
      # Convert GMIS values of 255 to missing
      map(., function(x) case_when(x == 200 ~ 0,
                                   x == 255 ~ NA_real_,
                                   TRUE ~ x)) %>% 
      map(., mean, na.rm = TRUE) %>%  # Ignore cells with missing GMIS values when calculating mean
      unlist()
    
    # Add column with GMIS values
    df_out <- df %>% 
      mutate(GMIS_MEAN = gmis_data_mod,
             GMIS_MEAN = ifelse(GMIS_MEAN == 255, NA, GMIS_MEAN))
    
    write_csv(df_out, path = outpath, col_names = TRUE)
  }
}

# Load all population-mean dataframes as list
inpath <- 'data/clean/popMeans_allCities/'
popMeans_dfList <- create_df_list(inpath)

# Map GMIS extraction function over dataframe list
inpath <- '~/Downloads/gmis_datasets/'
outpath <- 'data/raw/environmental_data/Extracted_GMIS_filtered/'
dir.create(outpath)
popMeans_subset <- popMeans_dfList[157]
purrr::walk(popMeans_dfList, extract_gmis, inpath = inpath, outpath = outpath)
