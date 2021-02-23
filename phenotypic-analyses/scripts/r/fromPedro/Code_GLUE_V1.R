# Script to perform robust regression for each city
#
# Author: Pedro Peres-Neto
# Sent by Email on March 27
rm(list=ls()) 

# Load required packages
library(MASS)
library(tidyverse)

### FUNCTION NOT CURRENTLY USED. WILL ANNOTATE LATER ###
dist.based.RDA <- function(D, X){
  D <- as.matrix(D)
  X <- as.matrix(X)
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
  
  # Projector matrix H
  X.c <- scale(X, center=TRUE, scale=FALSE)   # Centre matrix X
  m <- qr(X.c, tol=1e-6)$rank   
  # F statistic
  H <- (X.c[,1] %*% t(X.c[,1]))/((t(X.c[,1]) %*% X.c[,1])[1,1]) 
  HGH <- H %*% G %*% H
  SSYhat <- sum(diag(HGH))
  F <- SSYhat/(SSY-SSYhat)  
  Rsquare <- SSYhat/SSY     
  RsqAdj <- 1-((1-Rsquare)*(n-1)/(n-1-m))
  list(F=F*(n-m-1)/m, Rsquare=c(Rsquare,RsqAdj), PCoA.vectors=princ.coord)
}

#' Perform Principal Coordinate Analysis (PCoA)
#' 
#' @param D An M x N matrix
#' 
#' @return Scaled Eigenvalues from PCoA
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

#' Perform robust regression
#' 
#' @param data.reg Vector to be used as response variable
#' @param dist.UrbanCenter Vector to be used as predictor variable
#' 
#' @return List with three elements: (1) 'fitted_vals': Predicted values along
#'     each position of `dist.UrbanCenter`, (2) 'residual_vals': Residuals from 
#'     robust regression, (3) 'slope': Slope relating `data.reg` to `dist.UrbanCenter`.
run.rlm <- function(data.reg,dist.UrbanCenter){
  rlm.res <- rlm(data.reg~dist.UrbanCenter,na.action=na.omit)
  fitted_vals <- rlm.res$fitted.values
  residual_vals <- rlm.res$residuals
  slope <- rlm.res$coefficients[2]
  if(all(is.na(rlm.res$fitted.values))){
    lm.res <- lm(data.reg~dist.UrbanCenter,na.action=na.omit)
    fitted_vals <- lm.res$fitted.values
    residual_vals <- lm.res$residuals
    slope <- lm.res$coefficients[2]
  }
  result <- list(fitted_vals=fitted_vals,residual_vals=residual_vals,slope=slope)
  return(result)
}

#' Generate predicted values of environmental variable using robust regression
#' 
#' @param all.data Dataframe with population-mean HCN frequencies and environmental data,
#'     and population distances for every city.
#' @param permute Whether populations within a city should be randomly reordered for 
#'     permutation test
#' 
#' @return List with three elements: (1) 'UrbanNonUrbanExtremeSites': Dataframe with mean 
#'     value of each environmental variable in the most urban and rural sites. (2) 'all.pred':
#'     Dataframe with predicted value of environmental variable from robust regressions for
#'     each position along 'std_distance'. (3) 'slope.matrix': Dataframe with slopes from 
#'     robust regression of environmental variable against 'std_distance' for each variable 
#'     and city processed independently. 
generate.pred.values <- function(all.data,permute=FALSE){
  options(warn=-1)
  city.names <- unique(all.data[,1])
  n.cities <- length(city.names)
  all.pred <- c() # Dataframe to store results
  num_vars <- dim(all.data %>% dplyr::select(contains("Mean")))[2]
  UrbanNonUrbanExtremeSites <- matrix(0,2,num_vars) # 2 represents rural and urban
  UrbanNonUrbanExtremeSites.df <- c()
  slope.matrix <- c()
  tmp.matrix <- matrix(0,num_vars,1)
  # Iterate over cities
  for (i in 1:n.cities){
    # print(as.character(city.names)[i])
    rows.For.city <- as.integer(which(all.data[,1]==city.names[i]))
    # Select all environmental variables
    data.reg <- all.data %>% slice(rows.For.city) %>% dplyr::select(contains("Mean"))
    dist.UrbanCenter <- all.data[rows.For.city,"std_distance"] # most rural is 1
    n.sites <- length(dist.UrbanCenter)
    pred.city <- matrix(0, n.sites, num_vars)
    if (!permute){
      for (j in 1:num_vars){
        if(!all(is.na(data.reg[,j]))){
          res.run.rlm <- run.rlm(data.reg[,j],dist.UrbanCenter)
          fitted_vals <- res.run.rlm$fitted_vals
          fitted_vals <- replace(data.reg[,j], !is.na(data.reg[,j]), fitted_vals)
          tmp.matrix[j] <- res.run.rlm$slope
        }else{
          fitted_vals <- rep(NA, n.sites)
        }
        pred.city[,j] <- fitted_vals
      }
    } 
    if (permute){
      order.sites <- sample(n.sites)
      for (j in 1:num_vars){
        if(!all(is.na(data.reg[,j]))){
          res.run.rlm <- run.rlm(data.reg[,j],dist.UrbanCenter)
          y.rnd<-res.run.rlm$fitted_vals+res.run.rlm$residual_vals[order.sites]
          res.run.rlm <- run.rlm(y.rnd,dist.UrbanCenter)
          fitted_vals <- res.run.rlm$fitted_vals
          fitted_vals <- replace(data.reg[,j], !is.na(data.reg[,j]), fitted_vals)
        }else{
          fitted_vals <- rep(NA, n.sites)
        }
        pred.city[,j] <- fitted_vals
      }
    }
    tmp <- cbind(dist.UrbanCenter,pred.city)
    tmp <- na.omit(tmp[order(tmp[,1]),])[,-1]
    if(!all(is.na(tmp))){
      slope.matrix <- cbind(tmp.matrix,slope.matrix)
      UrbanNonUrbanExtremeSites[1,1:18] <- tmp[1,]
      UrbanNonUrbanExtremeSites[2,1:18] <- tmp[dim(tmp)[1],]
      colnames(pred.city) <- paste(names(data.reg), "rlmFitted", sep = "_")
      all.pred <- rbind(data.frame(pred.city, distance= dist.UrbanCenter, city=city.names[i]), all.pred)
      colnames(UrbanNonUrbanExtremeSites) <- paste(names(data.reg), "rlmFitted", sep = "_")
      UrbanNonUrbanExtremeSites.df <- rbind(data.frame(UrbanNonUrbanExtremeSites, zone=rbind("MostUrban","MostRural"), city=city.names[i]), UrbanNonUrbanExtremeSites.df)
    }
  } # end for loop for cities
  slope.matrix <- t(slope.matrix)
  options(warn=0)
  result <- list(UrbanNonUrbanExtremeSites=UrbanNonUrbanExtremeSites.df,all.pred=all.pred,slope.matrix=slope.matrix)
  return(result)
}

#' PCoA of Urban/Rural environmental values using 'Zone' (i.e. Urban or Rural) and City
#'     constraing (i.e. predictor) variables
#' 
#' @param all.data Dataframe with population-mean HCN frequencies and environmental data,
#'     and population distances for every city.
#' @param permute Whether Residuals of reduced model (i.e., no interaction between city and zone)
#'     should be resampled for permutation analysis. 
#' 
#' @return List with two or four elements, depending on `permute`: If `permute = FALSE`, function
#'     returns: (1)  'D.UR': Vector with distance between urban and rural habitats in multivariate environment
#'     space, (2) 'V.angle.UR': Scaled Eigenvalues from PCoA of pairwise multivariate environmental vector angles,
#'     (3) 'V.slopes': Scaled Eigenvalues from PCoA of robust regression slopes for each environmental variable,
#'     (4) 'stats': F-stats from MANOVA with urban/rural environmental matrix as Response and City, Urban/Rural as factors.
#'     If `permute = TRUE`, returns: (1) 'fitted.MANOVA' fitted values with from MANOVA with urban/rural environmental matrix as 
#'     Response and City, Urban/Rural as factor, (2) 'stats': F-stats from MANOVA with urban/rural environmental 
#'     matrix as Response and City, Urban/Rural as factors.
calculate.stats <- function(all.data,permute=FALSE){
  result.pred <- generate.pred.values(all.data,permute)
  slope.matrix <- result.pred$slope.matrix
  num_vars <- dim(all.data %>% dplyr::select(contains("Mean")))[2]
  y.mat <- scale(as.matrix(result.pred$UrbanNonUrbanExtremeSites[,1:num_vars])) # Response environmental vars in urban and rural pops
  # or PCoA <- P.Coord.A(dist(y.mat)); but y.mat doesn't have 0 eigenvalues so is fine. 
  # code cities and zone result.pred$UrbanNonUrbanExtremeSites[,19:20]
  colnames(y.mat) <- colnames(result.pred$UrbanNonUrbanExtremeSites[,1:num_vars])
  # sum(apply(y.mat, 2, function(x) any(is.na(x))))
  n.cities <- dim(y.mat)[1]/2
  
  # Encode city and urban/rural as factors
  x.mat <- as.matrix(cbind(as.factor(result.pred$UrbanNonUrbanExtremeSites[,(num_vars+1)]),
                           as.factor(result.pred$UrbanNonUrbanExtremeSites[,(num_vars+2)])))
  
  colnames(x.mat) <- c("zones","city")

  # add intercept and interaction terms:
  x.mat <- cbind(rep(1,dim(x.mat)[1]),x.mat,x.mat[,1]*x.mat[,2])
  
  if (permute){
    # reduced model (no interaction):
    fit <- manova(y.mat ~ x.mat[,2]+x.mat[,3])
    y.hat <-fit$fitted.values
    y.res<-y.mat-y.hat                            # Residuals of reduced model (to be permuted)
    y.rand<-y.hat+y.res[sample(dim(y.mat)[1]),]
    y.mat<-y.rand
  }
  
  # set up statistics; since each group (urban non-urban and cities are balanced in their observations)
  # we don't need to use least-square means, i.e., contrasts * slopes.full.model
  ### differences between urban and non-urban and cities environments (akin to phenotypic change vectors):
  # length of vectors:
  # 
  # MANOVA:
  fit <- manova(y.mat ~ x.mat[,2]*x.mat[,3]) # zone and city
  F.stats.manova <- summary(fit, test="Pillai")$stats[1:3,3] # F-stats
  pred.manova <- fit$fitted.values # RDA purposes
  
  ### differences between urban and non-urban environments for each city (akin to phenotypic change vectors):
  sites.MostUrban <- which(result.pred$UrbanNonUrbanExtremeSites$zone=="MostUrban")
  sites.MostRural <- which(result.pred$UrbanNonUrbanExtremeSites$zone=="MostRural")
  diff.vector <- as.matrix(y.mat[sites.MostUrban,1:num_vars]-y.mat[sites.MostRural,1:num_vars])
  
  # length of vectors (this considers only change betwee urban and non-urban within the same city):
  distance.vector <- sqrt(diag((diff.vector)%*%t(diff.vector)))

  # length of changes between urban and rural, summed across all cities for significance test (but PCoA can be performed)
  diff.lengths <- as.matrix(dist(distance.vector))
  sum.city.contrasts <- sum(diff.lengths)/2 # statistic to test if cities vary in the magnitude of change betwee urban and rural

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
  
  sum.city.angles <- sum(angle.matrix)/2 # statistic to test if cities vary in the direction of change betwee urban and rural
  # principal coordinate of angles and slopes
  if (!permute){
    PCoA.angle.matrix <- P.Coord.A(sqrt(angle.matrix))
    PCoA.slope.matrix <- P.Coord.A(dist(scale(slope.matrix)))
  }
  if (!permute){
    stats <- c(F.stats.manova,sum.city.contrasts,sum.city.angles)
    names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","MANOVA.F.interaction","magnitude.change","direction.change")
    result<-list(D.UR=distance.vector,V.angle.UR=PCoA.angle.matrix,V.slopes=PCoA.slope.matrix,stats=stats)}
  if (permute){
    stats <- c(F.stats.manova,sum.city.contrasts,sum.city.angles)
    names(stats) <- c("MANOVA.F.zone","MANOVA.F.city","MANOVA.F.interaction","magnitude.change","direction.change")
    result<-list(fitted.MANOVA=fit$fitted.values,stats=stats)}
  return(result)
}

### THIS CODE DOESN'T RUN ###
permutation.tests(all.data,n.perm=100){
  result.stats.obs <- calculate.stats(all.data,permute=FALSE)
  # test averages
  stats.obs <- result.stats.obs$stats
  stats.rnd <- matrix(0,n.perm,5)
  for (i in 1:n.perm){
    print(i)
    result.stats.rnd <- calculate.stats(all.data,permute=TRUE)
    stats.rnd.tmp <- result.stats.rnd$stats
    stats.rnd[i,] <- stats.rnd.tmp
  }
  pval <- (rowSums(apply(stats.rnd, 1, function(i) i >= stats.obs)) + 1)  / (n.perm + 1)
}

# Setting working directory
# Identify all CSV files
csv.files <- list.files(pattern="*.csv")
# Assign all population-mean datasets to global environment
# Combine all datasets together into "all.data"
# Swansea not included due to high numbers of NA
# Woodstock not included as distance from urban has NAs

all.data <- c()
for (i in 1:length(csv.files)){
  data <- assign(csv.files[i], read.csv(csv.files[i]))
  all.data <- rbind(all.data,data)
}
all.data <- data.frame(all.data)
  
# analyses:
result.stats.obs <- calculate.stats(all.data,permute=FALSE)
D.UR <- result.stats.obs$D.UR
V.angle.UR <- result.stats.obs$V.angle.UR
V.slopes <- result.stats.obs$V.slopes
fitted.MANOVA <- result.stats.obs$fitted.MANOVA # a PCA on this matrix is equivalent to an RDA
