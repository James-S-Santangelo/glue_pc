# Script for the GLUE project analyses
# Author: Pedro Peres-Neto
rm(list=ls())

# Load required packages
library(MASS)
library(tidyverse)
library(heplots)
library(candisc)
library(ggord) # install_github('fawda123/ggord')
library(InPosition)
library(factoextra)
library(FactoMineR)
library(ggpubr)

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
  rlm.res <- rlm(scale(data.reg)~scale(dist.UrbanCenter),na.action=na.omit)
  #rlm.res <- rlm(data.reg~dist.UrbanCenter,na.action=na.omit)
  
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
  num_vars <- dim(all.data %>% select(contains("Mean"),contains("freqHCN")))[2]
  slopes <- matrix(0,num_vars,1)
  Predicted.Values <- c() # Dataframe to store predicted values
  Original.Values <- c() # Dataframe to store original values that were used to predict values, i.e., passed conditions of NA, etc.
  slope.dataFrame <- c()
  # Loops over cities
  for (i in 1:n.cities){
    pick.city.rows <- as.integer(which(all.data[,"city"]==city.names[i]))
    # Select environmental variables
    data.reg <- all.data %>% slice(pick.city.rows) %>% select(contains("Mean"),contains("freqHCN"))
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
  num_vars <- dim(Predicted.Values %>% select(contains("Mean"),contains("freqHCN")))[2]
  var.names <- names(Predicted.Values %>% select(contains("Mean"),contains("freqHCN")))
  
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
    pred.tmp <- Predicted.Values %>% slice(pick.city.rows) %>% select(contains("Mean"),contains("freqHCN"))
    original.tmp <- Original.Values %>% slice(pick.city.rows) %>% select(contains("Mean"),contains("freqHCN"))
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
  slope.matrix <- slope.matrix %>% select(-contains("freqHCN"))
  
  Predicted.Values <- result.pred$Predicted.Values
  Predicted.Values <- Predicted.Values %>% select(-contains("freqHCN"))
  Original.Values <- result.pred$Original.Values
  Original.Values <- Original.Values %>% select(-contains("freqHCN"))
  
  # for MANOVA it needs to be two extreme values to test for interactions
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  UrbanPredExtremes <- as.matrix(ExtreValues$UrbanPredExtremes)
  RuralPredExtremes <- as.matrix(ExtreValues$RuralPredExtremes)
  city.names <- ExtreValues$city.names
  n.cities <- length(unique(city.names))
  var.names <- names(result.pred$Predicted.Values %>% select(contains("Mean")))
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
    F.stats.manova <- summary(fit, test="Pillai")$stats[1:3,3]
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
  if (number.extreme.sites > 1){
    stats.rnd <- matrix(0,n.perm,5)
  } else{stats.rnd <- matrix(0,n.perm,4)}
  for (i in 1:n.perm){
    print(i)
    result.stats.rnd <- calculate.stats(all.data,permute=TRUE,number.extreme.sites=number.extreme.sites)
    stats.rnd[i,] <- result.stats.rnd$stats
  }
  pval <- (rowSums(apply(stats.rnd, 1, function(i) i >= stats.obs)) + 1)  / (n.perm + 1)
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
  fviz_pca_biplot(pca.result,
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
  ggarrange(p1,p2)
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
  sd.scores.city <- data.frame(select(city.df, -starts_with("comp")) %>% group_by(city,area) %>% summarise_all(funs(sd)))
  n.rows <- nrow(sd.scores.city)
  sd.scores.city <- data.frame(sd.scores.city,eig=matrix(rep(unique(city.df  %>% select(starts_with("comp"))),each=n.rows),n.rows))
  
  cities <- unique(city.df$city)
  n.cities <- length(cities)
  n.dimensions <- ncol(city.df %>% select(contains("comp")))
  vector.diff <- matrix(0,1,n.dimensions)
  vector.cities <- matrix(0,1,n.cities)
  vector.chiSquare <- matrix(0,1,n.cities)
  for (i in 1:n.cities){
    pick.city <- which(sd.scores.city$city==cities[i])
    mult.eig <- matrix(as.numeric(as.matrix(select(sd.scores.city[pick.city,], starts_with("Dim"))))*
                         as.numeric(as.matrix(select(sd.scores.city[pick.city,], starts_with("eig")))),2,n.dimensions)
    vector.cities[i] <- sum(abs(mult.eig[1,]-mult.eig[2,]))
  }
  ##### correction via permutation
  vector.diff.perm <- matrix(0,1,n.dimensions)
  vector.cities.perm <- matrix(0,n.cities,n.sim)
  for (j in 1:n.sim){
    print(j)
    # no need to rerun the PCA because the standardization does not change the correlation structure of the dev. matrix
    city.df.perm <- perm.within.blocks(city.df)
    sd.scores.city.perm <- data.frame(select(city.df.perm, -starts_with("comp")) %>% group_by(city,area) %>% summarise_all(funs(sd)))
    sd.scores.city.perm <- data.frame(sd.scores.city.perm,eig=matrix(rep(unique(city.df  %>% select(starts_with("comp"))),each=n.rows),n.rows))
    for (i in 1:n.cities){
      pick.city <- which(sd.scores.city.perm$city==cities[i])
      mult.eig <- matrix(as.numeric(as.matrix(select(sd.scores.city.perm[pick.city,], starts_with("Dim"))))*
                           as.numeric(as.matrix(select(sd.scores.city.perm[pick.city,], starts_with("eig")))),2,n.dimensions)
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
  
  Predicted.Values <- result.pred$Predicted.Values %>% select(-contains("freqHCN"))
  Original.Values <- result.pred$Original.Values %>% select(-contains("freqHCN"))
  
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
  tmp <- combine.data.cent %>% select(contains("Mean"))
  colnames(tmp) <- paste(colnames(tmp),".ctr",sep = "")
  
  city.df.all <- data.frame(cbind(combine.data.df,area=combine.data.cent$area,city=combine.data.cent$city,tmp))
  
  return(city.df.all)
}  

multi.disp.analysis <- function(res.dist,plot.disp=FALSE){
  ### analyses
  # Box M test:
  box  <- boxM(as.matrix(select((res.dist %>% select(contains("Mean"))), -contains("ctr"))), factor(res.dist$area))
  # plot(box)
  
  # Levene's test
  data.dev.lm <- manova(as.matrix(res.dist %>% select(contains("ctr")))  ~ factor(res.dist$area))
  Anova(data.dev.lm)
  data.can <- candisc(data.dev.lm)
  plot(data.can, which=1)
  
  # PCA; to centre elipsoids use group.center(data.dev,groups):
  pca.result <- PCA(sqrt(as.matrix(res.dist %>% select(contains("ctr")))), graph = FALSE,scale.unit=TRUE)
  pca.result$groups <- factor(res.dist$area)
  if (plot.disp==TRUE){pca.plots.disp1(pca.result)
  pca.plots.disp2(pca.result)}
  
  # standardize variables according to their effect sizes, i.e., F statistic:
  std.dev <- std.var.size.effect(as.matrix(sqrt(as.matrix(res.dist %>% select(contains("ctr"))))),res.dist$area)
  n.vars <- ncol(as.matrix(select((res.dist %>% select(contains("Mean"))), -contains("ctr"))))
  pca.result.std.dev <- PCA(std.dev, graph = FALSE,scale.unit=FALSE,ncp=n.vars)
  pca.result.std.dev$groups <- factor(res.dist$area)
  if (plot.disp==TRUE){pca.plots.disp2(pca.result.std.dev)}
  
  # graph based on mean values of dispersion per area per city
  std.dev.mean.per.city <- aggregate(std.dev,by=list(area=res.dist$area,city=res.dist[,"city"]),FUN=mean)
  pca.result.std.dev <- PCA(std.dev.mean.per.city[,3:(n.vars+2)], graph = FALSE,scale.unit=FALSE,ncp=n.vars)
  pca.result.std.dev$groups <- factor(std.dev.mean.per.city[,"area"])
  if (plot.disp==TRUE){pca.plots.disp2(pca.result.std.dev)}
  
  # update city.df with pca.scores:
  res.dist <- cbind(res.dist,pca.result.std.dev$ind$coord,t(pca.result.std.dev$eig[,1]))
  city.contribution.res <- city.contribution.variation(res.dist,n.sim=100)
  
  if (plot.disp==TRUE){ggplot(city.contribution.res$cont.corrected, aes(y=contribution, x=city)) +
    geom_bar(position="dodge", stat="identity") +
    theme(axis.text.x = element_text(angle = 90),axis.text=element_text(size=5,face="bold"))}
  
  return(city.contribution=city.contribution.res$cont.corrected)
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

pca.plots <- function(all.data,number.extreme.sites=1){
  result.pred <- generate.pred.values(all.data,permute=FALSE)
  
  Predicted.Values <- result.pred$Predicted.Values %>% select(-contains("freqHCN"))
  Original.Values <- result.pred$Original.Values %>% select(-contains("freqHCN"))
  
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  UrbanPredExtremes <- as.matrix(ExtreValues$UrbanPredExtremes)
  RuralPredExtremes <- as.matrix(ExtreValues$RuralPredExtremes)
  UrbanOriginalExtremes <- as.matrix(ExtreValues$UrbanExtremes)
  RuralOriginalExtremes <- as.matrix(ExtreValues$RuralExtremes)
  
  # PCA based on the predicted values
  city.names <- unique(ExtreValues$city.names)
  n.cities <- length(unique(city.names))
  pca.res <- princomp(scale(rbind(UrbanPredExtremes,RuralPredExtremes)))
  Urban.pca <- pca.res$scores[1:n.cities,1:2]
  Rural.pca <- pca.res$scores[(n.cities+1):(2*n.cities),1:2]
  par(mfrow=c(1,2))
  parallel.line.plot(Urban.pca,Rural.pca)
  
  # PCA based on the original values
  pca.res <- princomp(scale(rbind(UrbanOriginalExtremes,RuralOriginalExtremes)))
  Urban.pca <- pca.res$scores[1:n.cities,1:2]
  Rural.pca <- pca.res$scores[(n.cities+1):(2*n.cities),1:2]
  parallel.line.plot(Urban.pca,Rural.pca)
  
  PCA.result<-princomp(rbind(UrbanOriginalExtremes,RuralOriginalExtremes),cor=TRUE)
  pdf("PCApbiplot.pdf",width=12,height=6,paper='special')
  biplot(PCA.result,xlabs=c(city.names,city.names))
  dev.off()
  
  valuesUrban <- data.frame(UrbanPredExtremes)
  valuesRural <- data.frame(RuralPredExtremes)
  rownames(valuesUrban)=as.character(city.names)
  #rownames(valuesRural)=as.character(rep("",n.cities))
  valuesUrbanRural <- rbind(valuesUrban,valuesRural)
  res.pca <- princomp(valuesUrbanRural, cor = TRUE)
  pdf("PCApbiplot1.pdf",width=12,height=6,paper='special')
  fviz_eig(res.pca)
  dev.off()
  pdf("PCApbiplot2.pdf",width=12,height=6,paper='special')
  fviz_pca_ind(res.pca,label=c("ind",city.labels),col.ind = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
  
  fviz_pca_ind(res.pca,col.ind = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
  
  dev.off()
  pdf("PCApbiplot3.pdf",width=12,height=6,paper='special')
  fviz_pca_var(res.pca,col.var = "contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),repel = TRUE)
  dev.off()
  pdf("PCApbiplot4.pdf",width=12,height=6,paper='special')
  fviz_pca_biplot(res.pca, repel = TRUE,col.var = "#2E9FDF",col.ind = "#696969")
  dev.off()
  
  groups <- as.factor(c(rep("Urban",each=n.cities),rep("Rural",each=n.cities)))
  
  pdf("PCApbiplot5.pdf",width=12,height=6,paper='special')
  fviz_pca_ind(res.pca,col.ind = groups,palette = c("#00AFBB",  "#FC4E07","#E7B800"),addEllipses = TRUE, ellipse.type = "convex", legend.title = "Groups",repel = TRUE)
  dev.off()
}

slope.freqHCN <- function(all.data,number.extreme.sites=1){
  # to make is compatible with number of complete environmental variable cases
  result.pred <- generate.pred.values(all.data,permute=FALSE)
  Predicted.Values <- result.pred$Predicted.Values
  Original.Values <- result.pred$Original.Values
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  city.names <- unique(ExtreValues$city.names)
  n.cities <- length(city.names)
  slopes <- matrix(0,n.cities,2) # 1: regular OLS slope & 2: robust slope
  # Loops over cities
  is.sign <- matrix(0,n.cities,1)
  for (i in 1:n.cities){
    pick.city.rows <- as.integer(which(all.data[,"city"]==city.names[i]))
    lm.res <- lm(scale(all.data[pick.city.rows,"freqHCN"])~scale(all.data[pick.city.rows,"std_distance"]))
    slopes[i,1] <- coefficients(lm.res)[2]
    slopes[i,2] <- run.rlm(scale(all.data[pick.city.rows,"freqHCN"]),scale(all.data[pick.city.rows,"std_distance"]))$slope
    if (anova(lm.res)$"Pr(>F)"[1] < 0.05) {is.sign[i] <- 1}
  }
  colnames(slopes) <- c("ols","rlm")
  result <- list(city=city.names,slopes=slopes,is.sign=is.sign)
  return(result)
  
  data.frame(cbind(city.names,slopes))
}

extreme.freqHCN <- function(all.data,number.extreme.sites=1){
  result.pred <- generate.pred.values(all.data,permute=FALSE)
  Predicted.Values <- result.pred$Predicted.Values
  Original.Values <- result.pred$Original.Values
  
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=number.extreme.sites)
  UrbanPredExtremes <- ExtreValues$UrbanPredExtremes[,"freqHCN"]
  RuralPredExtremes <- ExtreValues$RuralPredExtremes[,"freqHCN"]
  UrbanOriginalExtremes <- ExtreValues$UrbanExtremes[,"freqHCN"]
  RuralOriginalExtremes <- ExtreValues$RuralExtremes[,"freqHCN"]
  
  result <- list(UrbanOriginalExtremes=UrbanOriginalExtremes,RuralOriginalExtremes=RuralOriginalExtremes,pred.diff=UrbanPredExtremes-RuralPredExtremes,obs.diff=UrbanOriginalExtremes-RuralOriginalExtremes)
  return(result)
}

perm.std_dist.within.cities <- function(all.data){
  cities <- unique(all.data$city)
  n.cities <- length(cities)
  std_dist.perm <- all.data
  for (i in 1:n.cities){
    pick.city <- which(all.data$city==cities[i])
    std_dist.perm[pick.city,"std_distance"] <- sample(all.data[pick.city,"std_distance"])
  }
  return(std_dist.perm)
}

ReactionNorm.line.plot <- function(Urban,Rural){
  xlab="Urban -----> Rural"
  ylab="freqHCN"
  n <- length(Urban)
  X <- c(rep(1,n),rep(2,n))
  plot(X,c(Urban,Rural), xlab = xlab, ylab = ylab,
       pch = 21, bg = c(rep("red",n),rep("green",n)), col = "black",
       lwd = 1, cex = 0.75,xaxt="n")
  Y2 <- matrix(0,n,2)
  X2[1:n,1] <- X[1:n]
  X2[1:n,2] <- X[(n+1):(2*n)]
  Y2[1:n,1] <- Urban
  Y2[1:n,2] <- Rural
  for (i in 1:n){
    segments(X2[i,1],Y2[i,1],X2[i,2],Y2[i,2])
  }
}

avg.env.per.city <- function(all.data){
  # too keep the same cities used to general predictions and slopes
  result.pred <- generate.pred.values(all.data,permute=FALSE)
  Predicted.Values <- result.pred$Predicted.Values 
  Original.Values <- result.pred$Predicted.Values 
  ExtreValues <- pick.extreme.values(Predicted.Values,Original.Values,number.extreme.sites=1)
  city.names <- ExtreValues$city.names
  sub.set.data <- filter(all.data, city %in% city.names)
  var.names <- names(sub.set.data %>% select(contains("Mean")))
  mean.per.group <- sub.set.data %>% group_by(city) %>%
    summarise_at(vars(var.names), list(mean),na.rm=TRUE)
  mean.per.group <- data.frame(mean.per.group)
  return(mean.per.group)
}


# Set working directory
# Identify all CSV files
csv.files <- list.files(pattern="*.csv")
# Assign all population-mean datasets to global environment
# Combine all datasets together into "all.data"
# Swansea not included due to high numbers of NA
# Woodstock not included as distance from urban has NAs

all.data <- c()
for (i in 1:length(csv.files)){
  data <- assign(csv.files[i], read.csv(csv.files[i]) %>% dplyr::select(city, std_distance, freqHCN, matches("*Mean$")))
  all.data <- rbind(all.data,data)
}

####################
# HCN clines as James figures
city.names <- unique(all.data$city)
n.cities <- length(city.names)
pick.city.rows <- as.integer(which(all.data[,"city"]==city.names[1]))
plot(all.data$freqHCN[pick.city.rows]~all.data$std_distance[pick.city.rows],cex=0,xlab="Standardized Distance",ylab="Frequency of HCN",xlim=c(-0.1,1.1),ylim=c(-0.1,1))
sign.cases <- matrix(0,n.cities,1)
for (i in 1:n.cities){
  pick.city.rows <- as.integer(which(all.data[,"city"]==city.names[i]))
  lm.res <-  rlm(all.data$freqHCN[pick.city.rows]~all.data$std_distance[pick.city.rows],na.action=na.omit)
  if(all(is.na(lm.res$fitted.values))){
    lm.res <- lm(all.data$freqHCN[pick.city.rows]~all.data$std_distance[pick.city.rows],na.action=na.omit)
  }
  if (coefficients(lm.res)[2] < 0){abline(lm.res,col="red")}
  if (coefficients(lm.res)[2] > 0){abline(lm.res,col="blue")}
}

# analyses:
result.stats.obs <- calculate.stats(all.data,permute=FALSE,number.extreme.sites=2)
D.UR <- result.stats.obs$D.UR
V.angle.UR <- result.stats.obs$V.angle.UR # principal coordinate of angles
# no need to use the principal coordinates of slopes (V.slopes); good to test angles but not to model HCN
# V.slopes <- result.stats.obs$V.slopes # principal coordinate of slopes of environmental variation
env.slopes <- as.matrix(result.stats.obs$slope.matrix[,1:9])
# permutation.tests(all.data,n.perm=99,number.extreme.sites=2); all significant at 0.001
res.dist <- mult.dispersion(all.data,number.extreme.sites=5)
res.dist.metric <- multi.disp.analysis(res.dist,plot.disp=FALSE)
avg.env.res <- as.matrix(avg.env.per.city(all.data)[,-1]) # to remove the city names

# modelling slopes of HCN
var.names <- c("AI","PET","DEM","GMIS","NDSI","sLST","sNDVI","sLST","sNDVI")
colnames(env.slopes) <- var.names
slopes.HCN <- slope.freqHCN(all.data,number.extreme.sites=1)
var.names <- c("mAI","mPET","mDEM","mGMIS","mNDSI","msLST","msNDVI","msLST","msNDVI")
colnames(avg.env.res) <- var.names
# the second columns of slopes.HCN$slopes are the rlm slopes
df <- data.frame(HCNslopes=slopes.HCN$slopes[,2],env.slopes,avg.env.res)

# contributions of each set of predictors
summary(lm(HCNslopes ~ res.dist.metric$contribution,data=df)) # Adjusted R-squared:  0.002875 
summary(lm(HCNslopes ~ D.UR,data=df)) # Adjusted R-squared:  0.0273; p-value: 0.0267
summary(lm(HCNslopes ~ env.slopes,data=df)) #Adjusted R-squared:  0.05374; p-value: 0.05667
summary(lm(HCNslopes ~ avg.env.res,data=df)) #Adjusted R-squared:  0.06653; p-value: 0.03085

# contributions combined of all predictors
preds <- cbind(D.UR,env.slopes,avg.env.res)
df <- data.frame(HCNslopes=slopes.HCN$slopes[,2],preds) 
summary(lm(slopes.HCN$slopes[,2] ~ preds,data=df)) #Adjusted R-squared:  0.08356; p-value: 0.04706

# a model considering interactions among predictors
interactions.terms <- model.matrix( ~.^2, data=data.frame(scale(preds)))[,-1] # remove intercept

# LASSO model selection
library(glmnet)
cvfit <- cv.glmnet(interactions.terms,slopes.HCN$slopes[,2])
coef(cvfit, s = "lambda.min")
non.zero.pred <- which(coef(cvfit, s = "lambda.min")!=0)[2:length(which(coef(cvfit, s = "lambda.min")!=0))]-1
colnames(interactions.terms[,non.zero.pred])
# select non-zero slopes:
int <- interactions.terms[,non.zero.pred]
lm.inter <- lm(scale(slopes.HCN$slopes[,2]) ~ int)
summary(lm.inter); # Adjusted R-squared:  0.272

##### demonstrates that the LASSO is not inflating model selection
# permute distances within cities
perm.dist <- perm.std_dist.within.cities(all.data)
# recalculate HCN slopes
slopes.HCN.perm <- slope.freqHCN(perm.dist,number.extreme.sites=1)
vfit <- cv.glmnet(interactions.terms,slopes.HCN.perm$slopes[,2])
non.zero.pred <- which(coef(cvfit, s = "lambda.min")!=0)[2:length(which(coef(cvfit, s = "lambda.min")!=0))]-1
colnames(interactions.terms[,non.zero.pred])
int <- interactions.terms[,non.zero.pred]
lm.inter <- lm(scale(slopes.HCN.perm$slopes[,2]) ~ int)
summary(lm.inter)$r.squared

############# graphs of HCN frequency across cities
ext.valuesHCN <- extreme.freqHCN(all.data,number.extreme.sites=1)
Urban <- ext.valuesHCN$UrbanOriginalExtremes
Rural <- ext.valuesHCN$RuralOriginalExtremes
ReactionNorm.line.plot(Urban,Rural)

############# miscelania (for now):
# PCA on predicted extreme environmental values
pca.plots(all.data,number.extreme.sites=1)

pdf("PCAplot.pdf",width=12,height=6,paper='special')
pca.plot(all.data,number.extreme.sites=1)
dev.off()

############# permutation of std_dist

per.std_dist <- perm.std_dist.within.cities(all.data)
result.stats.obs <- calculate.stats(per.std_dist,permute=FALSE,number.extreme.sites=5)
V.angle.UR <- result.stats.obs$V.angle.UR # principal coordinate of angles
ext.HCN <- extreme.freqHCN(per.std_dist,number.extreme.sites=1)
summary(lm(log(abs(ext.HCN$pred.diff)) ~ V.angle.UR[,1:100]))
slopes.HCN <- slope.freqHCN(per.std_dist)
summary(lm(slopes.HCN$slopes[,1] ~ V.angle.UR[,1:100]))

# forward selection
library(olsrr)
data.df <- data.frame(y=sqrt(abs(slopes.HCN$slopes[,1])),pred=V.angle.UR[,1:50])
model <- lm(y~.,data=data.df)
StepF.res <- ols_step_forward_p(model,penter=0.05)
StepF.res
preds <- as.numeric(gsub("pred.", "", StepF.res$predictors))











