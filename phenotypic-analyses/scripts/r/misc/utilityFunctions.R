# This script contains functions that are called by other R scripts
# as part of the GLUE analysis pipeline


#' Convert degrees to radians
#' 
#' @param deg Angle in degrees
#' 
#' @return Angle in radians
deg2rad <- function(deg){
  return(deg*pi/180)
}


#' Calculates the geodesic distance between two points specified 
#'     by radian latitude/longitude using the Haversine formula (hf)
#'     
#' @param long1 Longitude of first point in decimal degrees
#' @param lat1 Latitude of first point in decimal degrees
#' @param long2 Longitude of second point in decimal degrees
#' @param lat2 Latitude of first point in decimal degrees
#' 
#' @return Distance between two points in kilometers
haversine <- function(long1, lat1, long2, lat2) {
  
  # Ensure Lats and Longs are in radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  
  # Calculate geodesic distance based on havesine formala
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  d = R * c
  return(d) # Distance in km
}

#' Standardizes variable between 0 and 1. 
#'
#' @param df Dataframe containing variable to be standardized as column
#' @param var Variable to standardize as character string
#' 
#' @return Dataframe with standardized variable. Name prefixed with 'std_'
std_var_zero_one <- function(df, var){
  
  varname <- paste0('std_', var)
  
  df_out <- df %>% 
    mutate(max = max(!!sym(var), na.rm = TRUE),
           min = min(!!sym(var), na.rm = TRUE),
           # min_dist = min(distance),
           var_out = ((!!sym(var) - min)/ (max - min)),
           !!varname := round(var_out, 4)) %>% 
    dplyr::select(-min, -max, -var_out)
  
  return(df_out)
}

#' Creates list with dataframes as elements for all CSVs in inpath
#' 
#' @param inpath Path to directory containing CSVs to load
#' 
#' @return df_list. List of dataframes, each from a different CSV
create_df_list <- function(inpath){
  
  # Get all csv files in inpath
  files <- dir(inpath, pattern = "*.csv")
  
  # read in all the files, appending the path before the filename
  df_list <- files %>%
    map(~ read_csv(file.path(inpath, .)))
  
  return(df_list)
}


#' Identifies continent from latitude and ongitude
#' 
#' https://stackoverflow.com/questions/21708488/get-country-and-continent-from-longitude-and-latitude-point-in-r
#' 
#' @param points Dataframe with longitude as first column, latitude as second
#' 
#' @return Vector with continent names
coords2continent = function(points){  
  
  # countriesSP <- getMap(resolution='low')
  countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  #indices$continent   # returns the continent (6 continent model)
  indices$REGION   # returns the continent (7 continent model)
  #indices$ADMIN  #returns country name
  #indices$ISO3 # returns the ISO3 code 
}

#' Filters raw environmental data into clean data for later merging
#'     with clean population-mean datasets
#'     
#' @param df Dataframe with raw environmental data to be filtered
#' @param variable Name ('character') of environmental variable.
#'     Used to rename variable in dataset.
#' @param outpath Path ('character') to which filtered environmental 
#'     data will be written.
#'     
#' @return None. Writes filtered data to disk in outpath. 
filter_environmental_data <- function(df, variable, outpath){
  
  # Get name of city being processed
  city_name <- df$city[1]

  # Edge case for Newhaven.
  if(city_name == "New_Haven"){
    city_name <- "Newhaven"
    df <- df %>% mutate(city = fct_recode(city, 'Newhaven' = 'New_Haven'))
  }
  
  # Retrieve last 2 columns (environmental data)
  enviro_data_filtered <- df %>%
    dplyr::select(city, population, tail(names(.), 1)) %>%
    
    # Rename variables
    rename_at(vars(contains('MEAN')), function(x) paste0(variable, "_", "Mean")) %>%
    
    # Forward fill empty cells
    fill(tail(names(.), 1), .direction = "down")
  
  # Full path to which data frame will be written
  path <- paste0(outpath, city_name, "_", variable, ".csv")
  # print(path)
  
  if(file.exists(outpath)){
    print(sprintf('%s data already exists for %s', variable, city_name))
  }
  
  # Write dataframe
  write_csv(enviro_data_filtered, path = path, col_names = TRUE)
}


#' Generates population-mean dataset from individual-plant data.
#' 
#' @param df Individual-plant dataframe for GLUE city
#' @param city_centers Dataframe with city center coodinates for
#'     calculating distances between points
#' @param outpath **optional** Path ('character') to which pop-mean dataset 
#'     should be written
#' 
#' @return None. Writes population-mean dataset to outpath
generate_popMeans <- function(df, city_centers, outpath){
  
  # Get name of city being processed
  city_name <- df$city[1]
  # print(city_name)
  
  # Create population-mean dataset
  df_popMeans <- df %>%
    mutate(hcn_result = as.numeric(hcn_result)) %>%
    
    # Group by city, lat/long to ensure columns are kept
    group_by(city, population, population_latitude, population_longitude) %>%
    
    # Calculate HCN frequency
    summarise(total_plants = sum(!is.na(hcn_result)),
              numCyanogenic = sum(hcn_result, na.rm = TRUE),
              freqHCN = round(numCyanogenic / total_plants, 3)) %>%
    filter(!is.nan(freqHCN)) %>%
    
    # Add city centers to pop-mean dataset
    left_join(., city_centers, by = "city") %>%
    
    # Use haversine to calculate distance to center
    mutate(distance = haversine(population_longitude, population_latitude,
                                longitude_city, latitude_city),
           distance = round(distance, 4)) %>%
    
    # Group by city to make sure column-wide max and min distance are retrieved
    group_by(city) %>%
    
    # Calculate standardized distance to center
    mutate(max_dist = max(distance),
           # min_dist = min(distance),
           std_distance = (distance / max_dist),
           std_distance = round(std_distance, 4))
    
    # Return dataframe if not writting to disk
    if(missing(outpath)){
      
      return(df_popMeans)
      
    }else{
      
      # Full path to which data frame will be written
      path <- paste0(outpath, city_name, ".csv")
      
      # Write dataframe
      write_csv(df_popMeans, path = path, col_names = TRUE)
  }
}

#' Clean MTJJ dataframes and standardize to GLUE dataset format
#' 
#' @param df Dataframe to be cleaned and processed
#' @param city_centers Dataframe with city center coodinates for
#'     calculating distances between points
#' @param outpath **optional** Path ('character') to which pop-mean dataset 
#'     should be written
#'     
#' @return df_popMeans. Cleaned and processed dataframe
clean_mtjj_popMeans <- function(df, city_centers, plantsPerPop, outpath){
  
  # Get name of city being processed
  city_name <- df$city[1]
  
  df_popMeans <- df %>%
    left_join(., city_centers, by = "city") %>%
    mutate(distance = haversine(population_longitude, population_latitude,
                                longitude_city, latitude_city)) %>%
    
    # Group by city to make sure column-wide max and min distance are retrieved
    group_by(city) %>%
    
    # Calculate standardized distance to center
    mutate(max_dist = max(distance, na.rm = TRUE),
           # min_dist = min(distance, na.rm = TRUE),
           std_distance = (distance / max_dist), 
           std_distance = round(std_distance, 4)) %>%
    rename(., "freqHCN" = "HCN_Result") %>% 
    
    # Add number of plants per pop
    left_join(., plantsPerPop, by = c('city', 'population')) %>% 
    
    # Select columns, preserving order of other pop mean datasets
    dplyr::select(city, population, population_latitude, population_longitude,
           total_plants, numCyanogenic, freqHCN, latitude_city, longitude_city,
           continent, distance, max_dist, std_distance, Country)
    
    
  # Return dataframe if not writting to disk
  if(missing(outpath)){

    return(df_popMeans)
    
  }else{
    
    # Full path to which data frame will be written
    path <- paste0(outpath, city_name, ".csv")

    # Write dataframe
    write_csv(df_popMeans, path = path, col_names = TRUE)
  }
}


#' Add environmental data columns to all clean and processed GLUE
#'     and MTJJ datasets
#'     
#' @param df Dataframe to which environmental data will be added
#' @param outpath Path ('character') to which final dataframe will be written
#' @param err_file Path to file to which errors will be written. Errors
#'     are the environmental variables that are missing for certain cities.
#'     
#' @return None. Writes final population-mean dataset with environmental
#'     data to disk.
add_enviro_data <- function(df, outpath, err_file){
  
  df <- df %>% 
    ungroup() %>% 
    mutate(population = as.character(population))
  
  # Directory with clean environmental data
  inpath = "data/clean/environmental_data"
  
  # Get name of city being processed
  city_name <- df$city[1]

  # Edge case for Newhave. Need to change name to load datasets
  if(city_name == "New_Haven"){
    city_name <- "Newhaven"
    df <- df %>% 
      mutate(city = "Newhaven")
  }
  
  # Load all environmental data for city into data frame list
  enviro_datasets <- dir(inpath, recursive = TRUE, full.names = TRUE, 
                         pattern = paste0(city_name, "_"))
  df_list <- lapply(enviro_datasets, read_csv, col_types = cols())
  
  
  cols <- c("annualAI_Mean", "annualPET_Mean", "DEM_Mean", "GMIS_Mean", "summerLST_Mean", 
            "summerNDVI_Mean", "winterLST_Mean", "winterNDVI_Mean", "NDSI_Mean")
  
  # Merge environmental data for city into single dataframe
  # Keep only columns with "Mean", which represent the actual values
  if(length(df_list) != 0){
    merged_enviro_data <- Reduce(function(...) left_join(...,
                                                         by = c("city", "population"),
                                                         all.x = TRUE), df_list) %>%
      # dplyr::select(city, population, ends_with("Mean")) %>% 
      mutate(population = as.character(population))
    
    missing_cols <- setdiff(cols, names(merged_enviro_data)) 

    merged_enviro_data[missing_cols] <- NA
    
    if(length(missing_cols) > 0){
      write_out <- paste0(length(cols) - length(missing_cols), 
                          " of 9 environmental datasets found for ", 
                          city_name, ". Missing:")
      write(write_out, file = err_file, append = TRUE)
      write(paste0("    ", missing_cols), file = err_file, append = TRUE)
    }
    
    popMean_df <- merged_enviro_data %>% 
      left_join(., df, by = c("city", "population")) %>%
      mutate_at(vars(ends_with("Mean")), 
                funs(scaled = scale(., center = TRUE, scale = TRUE))) %>% 
      dplyr::select(noquote(order(colnames(.)))) %>%
      dplyr::select(city, Country, continent, everything()) %>%
      rename("country" = "Country")
    
    # print(head(merged_enviro_data))
  }else{
    write_out <- paste0("0 of 9 environmental datasets found for ",
                        city_name, ". Missing:")
    write(write_out, file = err_file, append = TRUE)
    popMean_df <- df
    missing_cols <- setdiff(cols, names(popMean_df)) 
    popMean_df[missing_cols] <- NA
    popMean_df <- popMean_df %>%
      mutate_at(vars(ends_with("Mean")), 
                funs(scaled = scale(., center = TRUE, scale = TRUE))) %>% 
      dplyr::select(noquote(order(colnames(.)))) %>%
      dplyr::select(city, Country, continent, everything()) %>%
      rename("country" = "Country")
  }
  
  # print(head(popMean_df))
  # Full path to which data frame will be written
  path <- paste0(outpath, city_name, ".csv")
  # print(path)
  
  # Write dataframe
  write_csv(popMean_df, path = path, col_names = TRUE)
}


#' Calculates the mean of each environmental variable in dataframe
#' 
#' @param df Dataframe for which means will be calculated
#' 
#' @return Summary dataframe with means of environmental variables.
calculate_city_eviro_means <- function(df){
  
  # print(df$city[1])
  # print(missing_cols)
  enviro_summary <- df %>%
    group_by(city) %>%
    summarise(annualAI = mean(annualAI_Mean, na.rm = TRUE),
              annualPET = mean(annualPET_Mean, na.rm = TRUE),
              DEM = mean(DEM_Mean, na.rm = TRUE),
              GMIS = mean(GMIS_Mean, na.rm = TRUE),
              NDSI = mean(NDSI_Mean, na.rm = TRUE),
              summerLST = mean(summerLST_Mean, na.rm = TRUE),
              summerNDVI = mean(summerNDVI_Mean, na.rm = TRUE),
              winterLST = mean(winterLST_Mean, na.rm = TRUE),
              winterNDVI = mean(winterNDVI_Mean, na.rm = TRUE))
  
  return(enviro_summary)
} 

#' Calculate length of transect (in km)
#' 
#' @param df Population-mean HCN frequency dataset for city
#' 
#' @return Dataframe with 'city' and 'transect_length'
transect_length <- function(df){
  
  df_out <- df %>% 
    group_by(city) %>% 
    dplyr::select(city, distance) %>% 
    summarise(min_dist = round(min(distance, na.rm = TRUE), 3),
              max_dist = round(max(distance, na.rm = TRUE), 3),
              transect_length = round(max_dist - min_dist, 3))
  
  return(df_out)
  
}

#' Run first-order robust regression for each city in GLUE
#'   Response: population-mean HCN frequency
#'   Predictor: standardized distance to the urban core
#' 
#' @param df Population-mean HCN frequency dataset for city
#' 
#' @return Dataframe with summary of model output
rlmStats <- function(df, response){
  
  city <- df %>% distinct(city) %>% pull()
  response_var <- df %>% pull(response)
  
  # If response variable contains data
  if(!all(is.na(response_var))){
  
    # Run robust regression
    rlm_mod <- rlm(response_var ~ std_distance, data = df, maxit = 200)
    slope <- round(rlm_mod$coefficients["std_distance"], 3)

    # Robust F-tests from sfsmisc package (Wald test)
    ftest <- f.robftest(rlm_mod, var = "std_distance")
    pval <- round(ftest$p.value, 3)
  }else{
    slope <- NA
    pval <- NA
  }
  
  # Output dataframe
  df_out <- data.frame(city = city, betaRLM = slope, pvalRLM = pval, var = response)

  return(df_out)
}

#' Run logistic regression with HCN frequency as response and
#'   standardised distance as the sole predictor
#' 
#' @param df Population-mean HCN frequency dataset for city
#' 
#' @return Dataframe with summary of model output
logistic_regression_stats <- function(df, predictor){
  
  city <- df %>% distinct(city) %>% pull()
  pred <- df %>% pull(predictor)
  resp <- df %>% pull(freqHCN)
  
  # Perform logistic regression
  mod <- glm(resp ~ pred, weights = total_plants, data = df, family = 'binomial')
  
  betaLog <- round(summary(mod)$coefficients['pred', 'Estimate'], 3)
  pvalLog <- round(summary(mod)$coefficients['pred', 'Pr(>|z|)'], 3)
  yint <- round(summary(mod)$coefficients['(Intercept)', 'Estimate'], 3)
  
  df_out <- data.frame(city = city, yint = yint, betaLog = betaLog, pvalLog = pvalLog)
  
  return(df_out)
}


#' Calculate pariwise pearson correlations among columns in data matrix
#'     
#' @param Matrix Matrix containing variables of interest as columns
#'     
#' @return Dataframe with pairwise Pearson correlations in lower triangle
#'   and P-values in upper triangle
generate_correlation_df <- function(matrix){
  
  # Create correlation matrix
  corr_calc <- rcorr(matrix, type = "pearson")
  corrMat <- round(corr_calc$r, 3)
  
  # Assign P-values to upper triangle of correlation matrix
  corrMat[upper.tri(corrMat)] <- round(corr_calc$P[upper.tri(corr_calc$P)], 3)
  
  # Create dataframe
  corrMat <- as.data.frame(corrMat) %>%
    rownames_to_column()
  
  return(corrMat)
}

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}

# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
  abline(lm(y~x), col='red', lwd = 2)
}

#' Function to plot panels C, D, and E for figure 3
#' 
#' @param df Population-mean HCN frequency dataframe
#' 
#' @return ggplot object
fig2_hij <- function(df){
  
  city <- df %>% pull(city) %>% unique()
  col <- case_when(city == 'Muenster' ~ pal[3],
                   city == 'Freehold' ~ pal[2])
  
  plot <- ggplot(df, aes(x = std_distance, y = freqHCN, weight = total_plants)) +
    geom_point(size = 3.5) +
    geom_smooth(method = 'glm', 
                color = ifelse(city == "Temuco", "black", col),
                fill = ifelse(city == "Temuco", "grey", col),
                method.args = list(family = "binomial"),
                size = 1.5) +
    xlab("Standardized distance") +
    ylab("Frequency of HCN") +
    scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.1)) +
    coord_cartesian(ylim = c(-0.01, 1.01), xlim = c(0, 1)) +
    ng1
  
  return(plot)
}

#### FROM PEDRO ####

# Functions below were written by Pedro Peres-Neto for environmental analyses
# Function documentation my own

P.Coord.A <- function(D){
  D <- as.matrix(D)
  n <- nrow(D)
  # Gower centring:
  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  G <- -0.5 * mat %*% (D^2) %*% mat
  SSY <- sum(diag(G))
  # Principal coordinate analysis
  eig <- eigen(G, symmetric=TRUE)
  values <- eig$values     # All eigenvalues
  nonzero.PCoA.eig <- which(abs(values/SSY) > sqrt(.Machine$double.eps))
  values <- values[nonzero.PCoA.eig]
  vectors <- eig$vectors   # All eigenvectors, scaled to lengths 1
  select <- which(values > 0)
  princ.coord <- vectors[,select] %*% diag(sqrt(values[select]))
  return(princ.coord)
}

run.rlm <- function(data.reg,dist.UrbanCenter){
  options(warn=-1)
  # rlm.res <- rlm(scale(data.reg)~scale(dist.UrbanCenter),na.action=na.omit)
  rlm.res <- rlm(data.reg~dist.UrbanCenter,na.action=na.omit)
  
  fitted_vals <- rlm.res$fitted.values
  residual_vals <- rlm.res$residuals
  slope <- rlm.res$coefficients[2]
  if(all(is.na(rlm.res$fitted.values))){
    lm.res <- lm(data.reg~dist.UrbanCenter,na.action=na.omit)
    #lm.res <- lm(data.reg~dist.UrbanCenter,na.action=na.omit)
    fitted_vals <- lm.res$fitted.values
    residual_vals <- lm.res$residuals
    slope <- lm.res$coefficients[2]
  }
  options(warn=0)
  result <- list(fitted_vals=fitted_vals,residual_vals=residual_vals,slope=slope)
  return(result)
}

generate.pred.values <- function(all.data,permute=FALSE){
  city.names <- unique(all.data[,"city"])
  n.cities <- length(city.names)
  num_vars <- dim(all.data %>% dplyr::select(contains("Mean"),contains("freqHCN")))[2]
  slopes <- matrix(0,num_vars,1)
  Predicted.Values <- c() # Dataframe to store predicted values
  Original.Values <- c() # Dataframe to store original values that were used to predict values, i.e., passed conditions of NA, etc.
  slope.dataFrame <- c()
  # Loops over cities
  for (i in 1:n.cities){
    pick.city.rows <- as.integer(which(all.data[,"city"]==city.names[i]))
    # Select environmental variables
    data.reg <- all.data %>% slice(pick.city.rows) %>% dplyr::select(contains("Mean"),contains("freqHCN"))
    dist.UrbanCenter <- all.data[pick.city.rows,"std_distance"] # most rural is 1
    n.sites <- length(dist.UrbanCenter)
    pred.city <- matrix(0, n.sites, num_vars)
    colnames(pred.city) <- names(data.reg)
    if (!permute){
      for (j in 1:num_vars){
        if(!all(is.na(data.reg[,j]))){
          res.run.rlm <- run.rlm(data.reg[,j],dist.UrbanCenter) # res stands for result
          fitted_vals <- res.run.rlm$fitted_vals
          if (length(fitted_vals) < n.sites){
            comb.vectors <- cbind(data.reg[,j],dist.UrbanCenter)
            comb.vectors[!!rowSums(is.na(comb.vectors)),] <- NA
            comb.vectors <- comb.vectors[,1]
            fitted_vals <- replace(comb.vectors, !is.na(comb.vectors), fitted_vals)}
          slopes[j] <- res.run.rlm$slope
        }else{
          fitted_vals <- rep(NA, n.sites)
          slopes[j] <- NA
        }
        pred.city[,j] <- fitted_vals
      }
    }
    if (permute){
      for (j in 1:num_vars){
        if(!all(is.na(data.reg[,j]))){
          res.run.rlm <- run.rlm(data.reg[,j],dist.UrbanCenter) # res stands for result
          fitted_vals <- res.run.rlm$fitted_vals
          residual_vals <- res.run.rlm$residual_vals
          order.sites <- sample(length(residual_vals))
          y.rnd <- fitted_vals + residual_vals[order.sites]
          
          if (length(y.rnd) < n.sites){
            comb.vectors <- cbind(data.reg[,j],dist.UrbanCenter)
            comb.vectors[!!rowSums(is.na(comb.vectors)),] <- NA
            comb.vectors <- comb.vectors[,1]
            y.rnd <- replace(comb.vectors, !is.na(comb.vectors), y.rnd)}
          
          res.run.rlm <- run.rlm(y.rnd,dist.UrbanCenter)
          fitted_vals <- res.run.rlm$fitted_vals
          if (length(fitted_vals) < n.sites){
            comb.vectors <- cbind(data.reg[,j],dist.UrbanCenter)
            comb.vectors[!!rowSums(is.na(comb.vectors)),] <- NA
            comb.vectors <- comb.vectors[,1]
            fitted_vals <- replace(comb.vectors, !is.na(comb.vectors), fitted_vals)}
        }else{
          fitted_vals <- rep(NA, n.sites)
        }
        pred.city[,j] <- fitted_vals
      }
    } # end loop over environmental variables
    
    # for each city:
    if(!all(is.na(cbind(dist.UrbanCenter,pred.city)))){
      # slopes
      slope.dataFrame <- cbind(slope.dataFrame,slopes)
      # predicted data
      Predicted.Values <- rbind(Predicted.Values,data.frame(pred.city, std_distance = dist.UrbanCenter, city=city.names[i]))
      Original.Values <- rbind(Original.Values,data.frame(data.reg, std_distance = dist.UrbanCenter, city=city.names[i]))
    }
  } # end loop over cities
  slope.dataFrame <- t(slope.dataFrame)
  colnames(slope.dataFrame) <- names(data.reg)
  slope.dataFrame <- data.frame(slope.dataFrame, city=city.names)
  result <- list(Predicted.Values=Predicted.Values,Original.Values=Original.Values,slopes=slope.dataFrame)
  return(result)
}

pick.extreme.values <- function(Predicted.Values,Original.Values,number.extreme.sites=1){
  city.names <- unique(Predicted.Values[,"city"])
  n.cities <- length(city.names)
  num_vars <- dim(Predicted.Values %>% dplyr::select(contains("Mean"),contains("freqHCN")))[2]
  var.names <- names(Predicted.Values %>% dplyr::select(contains("Mean"),contains("freqHCN")))
  
  UrbanExtreme.predicted <- matrix(0,number.extreme.sites,num_vars)
  colnames(UrbanExtreme.predicted) <- var.names
  RuralExtreme.predicted <- matrix(0,number.extreme.sites,num_vars)
  colnames(RuralExtreme.predicted) <- var.names
  UrbanExtreme.original <- matrix(0,number.extreme.sites,num_vars)
  colnames(UrbanExtreme.original) <- var.names
  RuralExtreme.original <- matrix(0,number.extreme.sites,num_vars)
  colnames(RuralExtreme.original) <- var.names
  
  UrbanExtreme.predicted.dataFrame <- c()
  RuralExtreme.predicted.dataFrame <- c()
  UrbanExtreme.OriginalData.dataFrame <- c()
  RuralExtreme.OriginalData.dataFrame <- c()
  
  for (i in 1:n.cities){
    pick.city.rows <- as.integer(which(Predicted.Values[,"city"]==city.names[i]))
    pred.tmp <- Predicted.Values %>% slice(pick.city.rows) %>% dplyr::select(contains("Mean"),contains("freqHCN"))
    original.tmp <- Original.Values %>% slice(pick.city.rows) %>% dplyr::select(contains("Mean"),contains("freqHCN"))
    dist.UrbanCenter <- Predicted.Values[pick.city.rows,"std_distance"] # most rural is 1
    for (j in 1:num_vars){
      if(!all(is.na(pred.tmp[,j]))){
        tmp <- cbind(dist.UrbanCenter,pred.tmp[,j])
        tmp <- na.omit(tmp[order(tmp[,1]),])
        UrbanExtreme.predicted[1:number.extreme.sites,j] <- tmp[1:number.extreme.sites,2]
        tmp <- na.omit(tmp[order(tmp[,1],decreasing = TRUE),])
        RuralExtreme.predicted[1:number.extreme.sites,j] <- tmp[1:number.extreme.sites,2]
        
        tmp <- cbind(dist.UrbanCenter,original.tmp[,j])
        tmp <- na.omit(tmp[order(tmp[,1]),])
        UrbanExtreme.original[1:number.extreme.sites,j] <- tmp[1:number.extreme.sites,2]
        tmp <- na.omit(tmp[order(tmp[,1],decreasing = TRUE),])
        RuralExtreme.original[1:number.extreme.sites,j] <- tmp[1:number.extreme.sites,2]
      }else{
        UrbanExtreme.predicted[1:number.extreme.sites,j] <- rep(NA, number.extreme.sites)
        RuralExtreme.predicted[1:number.extreme.sites,j] <- rep(NA, number.extreme.sites)
        UrbanExtreme.original[1:number.extreme.sites,j] <- rep(NA, number.extreme.sites)
        RuralExtreme.original[1:number.extreme.sites,j] <- rep(NA, number.extreme.sites)
      }
    } # num_vars
    
    UrbanExtreme.predicted.dataFrame <- rbind(UrbanExtreme.predicted.dataFrame,(UrbanExtreme.predicted))
    RuralExtreme.predicted.dataFrame <- rbind(RuralExtreme.predicted.dataFrame,(RuralExtreme.predicted))
    UrbanExtreme.OriginalData.dataFrame <- rbind(UrbanExtreme.OriginalData.dataFrame,(UrbanExtreme.original))
    RuralExtreme.OriginalData.dataFrame <- rbind(RuralExtreme.OriginalData.dataFrame,(RuralExtreme.original))
  } # cities
  # keep cities that have values for all variables, i.e., no NA for a particular variable within a city
  keep.cities <- which(rowSums(is.na(UrbanExtreme.predicted.dataFrame)) == 0)
  UrbanExtreme.predicted.dataFrame <- as.matrix(UrbanExtreme.predicted.dataFrame[keep.cities,])
  RuralExtreme.predicted.dataFrame <- as.matrix(RuralExtreme.predicted.dataFrame[keep.cities,])
  UrbanExtreme.OriginalData.dataFrame <- as.matrix(UrbanExtreme.OriginalData.dataFrame[keep.cities,])
  RuralExtreme.OriginalData.dataFrame <- as.matrix(RuralExtreme.OriginalData.dataFrame[keep.cities,])

  city.names <- rep(city.names,each=number.extreme.sites)
  city.names <- city.names[keep.cities]

  result <- list(city.names=city.names,UrbanPredExtremes = UrbanExtreme.predicted.dataFrame,RuralPredExtremes = RuralExtreme.predicted.dataFrame,UrbanExtremes = UrbanExtreme.OriginalData.dataFrame,RuralExtremes = RuralExtreme.OriginalData.dataFrame)
  return(result)
}

calculate.stats <- function(all.data,permute=FALSE,number.extreme.sites=2){
  # make sure that package raster is not loaded (it generates issues here)
  result.pred <- generate.pred.values(all.data,permute)
  
  # keep cities that have values for all variables, i.e., no NA for a particular variable within a city
  slope.matrix <- result.pred$slopes[which(rowSums(is.na(result.pred$slopes)) == 0),]
  # remove freqHCN to just consider environmental variables
  slope.matrix <- slope.matrix %>% dplyr::select(-contains("freqHCN"))
  
  Predicted.Values <- result.pred$Predicted.Values
  Predicted.Values <- Predicted.Values %>% dplyr::select(-contains("freqHCN"))
  Original.Values <- result.pred$Original.Values
  Original.Values <- Original.Values %>% dplyr::select(-contains("freqHCN"))
  
  # for MANOVA it needs to be two extreme values to test for interactions
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  UrbanPredExtremes <- as.matrix(ExtreValues$UrbanPredExtremes)
  RuralPredExtremes <- as.matrix(ExtreValues$RuralPredExtremes)
  city.names <- ExtreValues$city.names
  n.cities <- length(unique(city.names))
  var.names <- names(result.pred$Predicted.Values %>% dplyr::select(contains("Mean")))
  data.MANOVA.pred <- data.frame(rbind(UrbanPredExtremes[,var.names],RuralPredExtremes[,var.names]),city=city.names,zone=c(rep("Urban",n.cities*number.extreme.sites),rep("Rural",n.cities*number.extreme.sites)))
  y.mat <- data.MANOVA.pred[,var.names]
  
  if (permute){
    # reduced model (no interaction):
    if (number.extreme.sites > 1){
      fit <- manova(as.matrix(y.mat) ~ zone+city,data.MANOVA.pred)
      # summary(fit, test="Pillai")
      y.rand <- fit$fitted.values+fit$residuals[sample(dim(y.mat)[1]),]
      y.mat <- y.rand} else{
        y.mat <- y.mat[sample(dim(y.mat)[1]),]
      }
  }
  # set up statistics; since each group (urban non-urban and cities are balanced in their observations)
  # we don't need to use least-square means, i.e., contrasts * slopes.full.model
  
  # MANOVA:
  if (number.extreme.sites > 1){
    fit <- manova(as.matrix(y.mat) ~ zone*city,data.MANOVA.pred) # zone and city
    #print(summary(fit, test="Pillai"))
    F.stats.manova <- summary(fit, test="Pillai")$stats[1:3,3]
    #print(F.stats.manova)
  } else{
    print("number of extreme sites equals 1; interaction can't be tested; MANOVA reduced to main effects")
    fit <- manova(as.matrix(y.mat) ~ zone+city,data.MANOVA.pred) # zone and city
    F.stats.manova <- summary(fit, test="Pillai")$stats[1:2,3]
  }
  # pred.manova <- fit$fitted.values # RDA purposes
  
  # calculate means of extremes
  # https://statisticsglobe.com/mean-by-group-in-r
  
  data.MANOVA.pred[,var.names] <- y.mat
  # it places zone and city in the first two columns so that variables are now from column 3 to 11
  mean.per.group <- data.MANOVA.pred %>% group_by(zone,city) %>%
    summarise_at(vars(var.names), list(mean))
  mean.per.group <- data.frame(mean.per.group)
  
  # diff.vector is akin to "phenotypic change vector"
  diff.vector <- as.matrix(filter(mean.per.group, zone == "Urban")[,var.names]-filter(mean.per.group, zone == "Rural")[,var.names])
  
  # length of vectors (this considers only change between urban and non-urban within the same city):
  distance.vector <- sqrt(diag((diff.vector)%*%t(diff.vector)))
  
  # distance.vector could be also constructed as:
  # tmp <- as.matrix(dist(rbind(y.mat[sites.MostUrban,],y.mat[sites.MostRural,])))
  # diag(tmp[1:144,145:288])
  # because we have 144 cities then tmp[1,145] = distance.vector[1], tmp[2,146] = distance.vector[2], and so on
  
  # length of changes between urban and rural, summed across all cities for significance test (but PCoA can be performed)
  diff.lengths <- as.matrix(dist(distance.vector))
  sum.city.contrasts <- sum(diff.lengths)/(n.cities*n.cities) # statistic to test if cities vary in the magnitude of change betwee urban and rural
  
  # direction of change - we can probably solve this using matrix multiplication (will work on that later on)
  angle.matrix <- matrix(0,n.cities,n.cities)
  for (i in 1:n.cities){
    for (j in 1:n.cities){
      if (j > i){
        angle.matrix[i,j] <- acos(t((diff.vector[i,])/as.vector(distance.vector[i]))%*%((diff.vector[j,])/as.vector(distance.vector[j])))
        angle.matrix[i,j] <- angle.matrix[i,j]*180/pi
        angle.matrix[j,i] <- angle.matrix[i,j]
      }
    }
  }
  sum.city.angles <- sum(angle.matrix)/(n.cities*n.cities) # statistic to test if cities vary in the direction of change betwee urban and rural
  if (!permute){
    # principal coordinate of angles and slopes
    
    PCoA.angle.matrix <- P.Coord.A(sqrt(angle.matrix))
    PCoA.slope.matrix <- P.Coord.A(dist(scale(slope.matrix[,var.names])))
    
    stats <- c(F.stats.manova,sum.city.contrasts,sum.city.angles)
    if (number.extreme.sites>=2){
      names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","MANOVA.F.interaction","magnitude.change","direction.change")
    } else {
      names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","magnitude.change","direction.change")
    }
    result<-list(slope.matrix=slope.matrix,D.UR=distance.vector,V.angle.UR=PCoA.angle.matrix,V.slopes=PCoA.slope.matrix,fitted.MANOVA=fit$fitted.values,stats=stats,city=unique(city.names))
  }
  if (permute){
    stats <- c(F.stats.manova,sum.city.contrasts,sum.city.angles)
    stats <- c(F.stats.manova,sum.city.contrasts,sum.city.angles)
    if (number.extreme.sites>=2){
      names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","MANOVA.F.interaction","magnitude.change","direction.change")
    } else {
      names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","magnitude.change","direction.change")
    }
    result<-list(stats=stats)}
  return(result)
}

permutation.tests <- function(all.data,n.perm=99,number.extreme.sites=2){
  result.stats.obs <- calculate.stats(all.data,permute=FALSE,number.extreme.sites=number.extreme.sites)
  # test averages
  stats.obs <- result.stats.obs$stats
  print("Observed F-stats of multivariate mean environmental change")
  print(stats.obs)
  if (number.extreme.sites > 1){
    stats.rnd <- matrix(0,n.perm,5)
  } else{stats.rnd <- matrix(0,n.perm,4)}
  for (i in 1:n.perm){
    print(i)
    result.stats.rnd <- calculate.stats(all.data,permute=TRUE,number.extreme.sites=number.extreme.sites)
    stats.rnd[i,] <- result.stats.rnd$stats
  }
  pval <- (rowSums(apply(stats.rnd, 1, function(i) i >= stats.obs)) + 1)  / (n.perm + 1)
  return(pval)
}

group.center <- function(data.mat,groups) {
  n.var <- ncol(data.mat)
  n <- nrow(data.mat)
  data.transformed <- matrix(0,n,n.var)
  for (i in 1:n.var){
    data.transformed[,i] <- data.mat[,i]-tapply(data.mat[,i],groups,mean,na.rm=T)[groups]
  }
  colnames(data.transformed) <- colnames(data.mat)
  return(data.transformed)
}

pca.plots.disp1 <- function(pca.result){
  plot <- fviz_pca_biplot(pca.result,
                  # Individuals
                  geom.ind = "point",
                  fill.ind = pca.result$groups, col.ind = "black",
                  pointshape = 21, pointsize = 2,
                  palette = "jco",
                  addEllipses = TRUE,
                  # Variables
                  alpha.var ="contrib", col.var = "contrib",
                  gradient.cols = "RdYlBu",
                  legend.title = list(fill = "groups", color = "Contrib",
                                      alpha = "Contrib")
  )+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  print(plot)
}

pca.plots.disp2 <- function(pca.result){
  p1<-fviz_pca_ind(pca.result,
                   geom.ind = "point", # show points only (nbut not "text")
                   col.ind = pca.result$groups, # color by groups
                   palette = "jco",
                   addEllipses = TRUE, # Concentration ellipses
                   legend.title = "groups")
  
  p2<-fviz_pca_var(pca.result, col.var = "contrib",
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
  )
  print(ggarrange(p1,p2))
}

std.var.size.effect <- function(data.mat,groups){
  n.var <- ncol(data.mat)
  F.values <- matrix(0,n.var,1)
  data.std.effect.size <- scale(data.mat)
  for (i in 1:n.var){
    data.dev.lm <- lm(data.mat[,i]  ~ factor(groups))
    F.values[i] <- Anova(data.dev.lm)[1,3]
    data.std.effect.size[,i] <- data.std.effect.size[,i]*F.values[i]
  }
  return(data.std.effect.size)
}

perm.within.blocks <- function(city.df){
  cities <- unique(city.df$city)
  n.cities <- length(cities)
  city.df.perm <- city.df
  for (i in 1:n.cities){
    pick.city <- which(city.df$city==cities[i])
    city.df.perm[pick.city,"area"] <- sample(city.df[pick.city,"area"])
  }
  return(city.df.perm)
}

city.contribution.variation <- function(city.df,n.sim=100){
  sd.scores.city <- data.frame(dplyr::select(city.df, -starts_with("comp")) %>% group_by(city,area) %>% summarise_all(funs(sd)))
  n.rows <- nrow(sd.scores.city)
  sd.scores.city <- data.frame(sd.scores.city,eig=matrix(rep(unique(city.df  %>% dplyr::select(starts_with("comp"))),each=n.rows),n.rows))
  
  cities <- unique(city.df$city)
  n.cities <- length(cities)
  n.dimensions <- ncol(city.df %>% dplyr::select(contains("comp")))
  vector.diff <- matrix(0,1,n.dimensions)
  vector.cities <- matrix(0,1,n.cities)
  vector.chiSquare <- matrix(0,1,n.cities)
  for (i in 1:n.cities){
    pick.city <- which(sd.scores.city$city==cities[i])
    mult.eig <- matrix(as.numeric(as.matrix(dplyr::select(sd.scores.city[pick.city,], starts_with("Dim"))))*
                         as.numeric(as.matrix(dplyr::select(sd.scores.city[pick.city,], starts_with("eig")))),2,n.dimensions)
    vector.cities[i] <- sum(abs(mult.eig[1,]-mult.eig[2,]))
  }
  ##### correction via permutation
  vector.diff.perm <- matrix(0,1,n.dimensions)
  vector.cities.perm <- matrix(0,n.cities,n.sim)
  for (j in 1:n.sim){
    print(j)
    # no need to rerun the PCA because the standardization does not change the correlation structure of the dev. matrix
    city.df.perm <- perm.within.blocks(city.df)
    sd.scores.city.perm <- data.frame(dplyr::select(city.df.perm, -starts_with("comp")) %>% group_by(city,area) %>% summarise_all(funs(sd)))
    sd.scores.city.perm <- data.frame(sd.scores.city.perm,eig=matrix(rep(unique(city.df  %>% dplyr::select(starts_with("comp"))),each=n.rows),n.rows))
    for (i in 1:n.cities){
      pick.city <- which(sd.scores.city.perm$city==cities[i])
      mult.eig <- matrix(as.numeric(as.matrix(dplyr::select(sd.scores.city.perm[pick.city,], starts_with("Dim"))))*
                           as.numeric(as.matrix(dplyr::select(sd.scores.city.perm[pick.city,], starts_with("eig")))),2,n.dimensions)
      vector.cities.perm[i,j] <- sum(abs(mult.eig[1,]-mult.eig[2,]))
    }
  }
  mean.perm <- apply(vector.cities.perm,1,mean)
  std.perm <- apply(vector.cities.perm,1,sd)
  vector.cities.corrected <- (vector.cities-mean.perm)/std.perm
  vector.cities <- data.frame(city=cities,contribution=t(vector.cities))
  vector.cities.corrected <- data.frame(city=cities,contribution=t(vector.cities.corrected))
  vector.chiSquare <- data.frame(city=cities,contribution=t(vector.chiSquare))
  return(list(cont=vector.cities,cont.corrected=vector.cities.corrected,chiSquare=vector.chiSquare))
}

mult.dispersion <- function(all.data,number.extreme.sites=5){
  result.pred <- generate.pred.values(all.data,permute=FALSE)
  
  Predicted.Values <- result.pred$Predicted.Values %>% dplyr::select(-contains("freqHCN"))
  Original.Values <- result.pred$Original.Values %>% dplyr::select(-contains("freqHCN"))
  
  ExtremeValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  UrbanPredExtremes <- data.frame(city.names=ExtremeValues$city.names,ExtremeValues$UrbanPredExtremes)
  RuralPredExtremes <- data.frame(city.names=ExtremeValues$city.names,ExtremeValues$RuralPredExtremes)
  UrbanOriginalExtremes <- data.frame(city.names=ExtremeValues$city.names,ExtremeValues$UrbanExtremes)
  RuralOriginalExtremes <- data.frame(city.names=ExtremeValues$city.names,ExtremeValues$RuralExtremes)
  
  city.names <- unique(ExtremeValues$city.names)
  n.cities <- length(unique(city.names))
  area <- c(rep("Urban",number.extreme.sites),rep("Rural",number.extreme.sites))
  # center with urban and rural for each city separately
  combine.data.cent <-c()
  combine.data.df <-c()
  original <- FALSE # if FALSE, it uses predicted values
  for (i in 1:n.cities){
    if (original==TRUE){
      pick.city.rows <- which(UrbanOriginalExtremes[,"city.names"]==city.names[i])
      combine.data <- rbind(UrbanOriginalExtremes[pick.city.rows,-1],RuralOriginalExtremes[pick.city.rows,-1])
      combine.data.cent <- rbind(combine.data.cent,data.frame(abs(group.center(as.matrix(combine.data),factor(area))),area=area,city=rep(city.names[i],number.extreme.sites*2)))
      combine.data.df <- rbind(combine.data.df,combine.data)
    }
    if (original==FALSE){
      pick.city.rows <- which(UrbanPredExtremes[,"city.names"]==city.names[i])
      combine.data <- rbind(UrbanPredExtremes[pick.city.rows,-1],RuralPredExtremes[pick.city.rows,-1])
      combine.data.cent <- rbind(combine.data.cent,data.frame(abs(group.center(as.matrix(combine.data),factor(area))),area=area,city=rep(city.names[i],number.extreme.sites*2)))
      combine.data.df <- rbind(combine.data.df,combine.data)
    }
  }
  tmp <- combine.data.cent %>% dplyr::select(contains("Mean"))
  colnames(tmp) <- paste(colnames(tmp),".ctr",sep = "")
  
  city.df.all <- data.frame(cbind(combine.data.df,area=combine.data.cent$area,city=combine.data.cent$city,tmp))
  
  return(city.df.all)
}  

multi.disp.analysis <- function(res.dist,plot.disp=FALSE){
  ### analyses
  # Box M test:
  box  <- boxM(as.matrix(dplyr::select((res.dist %>% dplyr::select(contains("Mean"))), -contains("ctr"))), factor(res.dist$area))
  # plot(box)
  
  # Levene's test
  data.dev.lm <- manova(as.matrix(res.dist %>% dplyr::select(contains("ctr")))  ~ factor(res.dist$area))
  Anova(data.dev.lm)
  data.can <- candisc(data.dev.lm)
  # plot(data.can, which=1)
  
  # PCA; to centre elipsoids use group.center(data.dev,groups):
  pca.result <- PCA(sqrt(as.matrix(res.dist %>% dplyr::select(contains("ctr")))), graph = FALSE,scale.unit=TRUE)
  pca.result$groups <- factor(res.dist$area)
  if (plot.disp==TRUE){pca.plots.disp1(pca.result)
    pca.plots.disp2(pca.result)}
  
  # standardize variables according to their effect sizes, i.e., F statistic:
  std.dev <- std.var.size.effect(as.matrix(sqrt(as.matrix(res.dist %>% dplyr::select(contains("ctr"))))),res.dist$area)
  n.vars <- ncol(as.matrix(dplyr::select((res.dist %>% dplyr::select(contains("Mean"))), -contains("ctr"))))
  pca.result.std.dev <- PCA(std.dev, graph = FALSE,scale.unit=FALSE,ncp=n.vars)
  pca.result.std.dev$groups <- factor(res.dist$area)
  if (plot.disp==TRUE){pca.plots.disp2(pca.result.std.dev)}
  
  # graph based on mean values of dispersion per area per city
  std.dev.mean.per.city <- aggregate(std.dev,by=list(area=res.dist$area,city=res.dist[,"city"]),FUN=mean)
  pca.result.std.dev <- PCA(std.dev.mean.per.city[,3:(n.vars+2)], graph = FALSE,scale.unit=FALSE,ncp=n.vars)
  pca.result.std.dev$groups <- factor(std.dev.mean.per.city[,"area"])
  pca.result.std.dev$city <- factor(std.dev.mean.per.city[,"city"])
  if (plot.disp==TRUE){pca.plots.disp2(pca.result.std.dev)}
  
  # update city.df with pca.scores:
  res.dist <- cbind(res.dist,pca.result.std.dev$ind$coord,t(pca.result.std.dev$eig[,1]))
  city.contribution.res <- city.contribution.variation(res.dist,n.sim=100)

  if (plot.disp==TRUE){ggplot(city.contribution.res$cont.corrected, aes(y=contribution, x=city)) +
      geom_bar(position="dodge", stat="identity") +
      theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=5,face="bold"))}
  
  return(list(boxM = box, levene_lm = data.dev.lm, levene_can = data.can, pca_stand = pca.result.std.dev, city.contribution=city.contribution.res$cont.corrected))
}

parallel.line.plot <- function(X1,X2,xlab="PC-1",ylab="PC-2"){
  xlim.c=c(min(c(X1[,1],X2[,1]))+min(c(X1[,1],X2[,1]))/10,max(c(X1[,1],X2[,1]))+max(c(X1[,1],X2[,1]))/10)
  ylim.c=c(min(c(X1[,2],X2[,2]))+min(c(X1[,2],X2[,2]))/10,max(c(X1[,2],X2[,2]))+max(c(X1[,2],X2[,2]))/10)
  plot(X1[,1],X1[,2],col="red",xlab=xlab,ylab=ylab,las = 1,pch=16,cex=0.75,xlim=xlim.c,ylim=ylim.c)
  points(X2[,1],X2[,2],col="green",pch=16,cex=0.75,xlim=xlim.c,ylim=ylim.c)
  for (i in 1:dim(X1)[1]){
    segments(X1[i,1],X1[i,2],X2[i,1],X2[i,2])
  }
}