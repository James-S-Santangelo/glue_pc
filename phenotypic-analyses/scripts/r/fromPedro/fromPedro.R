# Script to perform robust regression for each city
#
# Author: Pedro Peres-Neto
# Sent by Email on March 27

# Setting working directory not necessary here due to use of Rproject
# setwd("~/Dropbox/Manuscripts & Data - collaborations (not with lab people)/GLUE - Urban")

# Identify all CSV files
path = "data/clean/popMeans_allCities_withEnviro/"
csv.files <- list.files(path, pattern="*.csv")
# Swansea not included due to high numbers of NA
# Woodstock not included as distance from urban has NAs

# Load required packages
library(MASS)
library(tidyverse)

# Assign all population-mean datasets to global environment
# Combine all datasets together into "all.data"
all.data <- c()
for (i in 1:length(csv.files)){
  data <- assign(csv.files[i], read.csv(paste0(path, csv.files[i])))
  all.data <- rbind(all.data,data)
}


city.names <- unique(all.data[,1])
n.cities <- length(city.names)
all.pred <- c() # Dataframe to store results
# Iterate over cities
for (i in 1:n.cities){
  print(as.character(city.names)[i])
  rows.For.city <- which(all.data[,1]==city.names[i])
  # Select all environmental variables, both scaled and unscaled
  data.reg <- all.data %>% slice(rows.For.city) %>% select(contains("Mean"))
  # data.reg <- all.data[rows.For.city,4:9] # just picked the first 6 variables; discuss with James about which ones to use in the final analyses
  dist.UrbanCenter <- all.data[rows.For.city,"std_distance"]
  # dist.UrbanCenter <- (max(dist.UrbanCenter)-dist.UrbanCenter)/(max(dist.UrbanCenter)-min(dist.UrbanCenter))
  n.sites <- length(dist.UrbanCenter)
  num_vars <- length(data.reg)
  pred.city <- matrix(0, n.sites, num_vars)
  for (j in 1:num_vars){
    # print(j)
    # print(colname)
    if(!all(is.na(data.reg[,j]))){
      rlm.res <- rlm(data.reg[,j]~dist.UrbanCenter)
      fitted_vals <- rlm.res$fitted.values
      fitted_vals <- replace(data.reg[,j], !is.na(data.reg[,j]), fitted_vals)
    }else{
      fitted_vals <- rep(NA, n.sites)
    }
    pred.city[,j] <- fitted_vals
  }
  colnames(pred.city) <- paste(names(data.reg), "rlmFitted", sep = "_")
  all.pred <- rbind(data.frame(pred.city, distance= dist.UrbanCenter, city=city.names[i]), all.pred)
}

write_csv(all.pred, "data/clean/rlm_fittedValues.csv")
