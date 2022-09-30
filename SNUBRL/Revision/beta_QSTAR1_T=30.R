## Coefficients for BA and CNT (QSTAR, seasnality)

#rm(list=ls())
library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(caret)
library(patchwork)
library(splines2)
library(fields)


# modify your working directory with the data file if necessary (using the setwd() command):
#getwd()
load("data_train.RData")
data_odd_DF <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)
train_DF_full <- subset(data_odd_DF, year == 1995 & month == 5)
loc <- cbind(train_DF_full$lon, train_DF_full$lat)
n <- nrow(loc)




# weight matrix (cut-off value is 800km)
#d <- matrix(rep(0, n * n), nr = n, nc = n)
W_star <- Z <- matrix(rep(0, 9 * n^2), nr = 3 * n, nc = 3 * n)
S <- c()
scale <- 270    # (distance) scale 
T <- 30            
T_scale <- T/scale  # scaled time distance (2 year) 
for (i in 2:n){  
  for (j in 1:(i - 1)){
    dist <- distGeo(loc[i, ], loc[j, ])/1000  # distance between ith and jth locations (km)
    dist_scale <- dist/scale
    if (dist < sqrt(4800^2 + 4 * T^2)/6){   # 제일 먼 거리의 1/6을 cutoff로
      Z[i, j] <- exp(-dist_scale)
      Z[j, i] <- exp(-dist_scale)
      
      Z[i + n, j + n] <- exp(-dist_scale)
      Z[j + n, i + n] <- exp(-dist_scale)
      
      Z[i + 2 * n, j + 2 * n] <- exp(-dist_scale)
      Z[j + 2 * n, i + 2 * n] <- exp(-dist_scale)
      
      Z[i, j + n] <- exp(-sqrt(T_scale^2 + dist_scale^2)) 
      Z[j + n, i] <- exp(-sqrt(T_scale^2 + dist_scale^2)) 
      
      Z[j, i + n] <- exp(-sqrt(T_scale^2 + dist_scale^2))   
      Z[i + n, j] <- exp(-sqrt(T_scale^2 + dist_scale^2))   
      
      Z[i, j + 2 * n] <- exp(-sqrt(4 * T_scale^2 + dist_scale^2)) 
      Z[j + 2 * n, i] <- exp(-sqrt(4 * T_scale^2 + dist_scale^2)) 
      
      Z[j, i + 2 * n] <- exp(-sqrt(4 * T_scale^2 + dist_scale^2)) 
      Z[i + 2 * n, j] <- exp(-sqrt(4 * T_scale^2 + dist_scale^2)) 
      
      Z[i + n, j + 2 * n] <- exp(-sqrt(T_scale^2 + dist_scale^2))
      Z[j + 2 * n, i + n] <- exp(-sqrt(T_scale^2 + dist_scale^2))
      
      Z[j + n, i + 2 * n] <- exp(-sqrt(T_scale^2 + dist_scale^2))
      Z[i + 2 * n, j + n] <- exp(-sqrt(T_scale^2 + dist_scale^2))
    } else {
      Z[i, j] <- 0
      Z[j, i] <- 0
      
      Z[i + n, j + n] <- 0
      Z[j + n, i + n] <- 0
      
      Z[i + 2 * n, j + 2 * n] <- 0
      Z[j + 2 * n, i + 2 * n] <- 0
      
      Z[i, j + n] <- 0  
      Z[j + n, i] <- 0 
      
      Z[j, i + n] <- 0 
      Z[i + n, j] <- 0
      
      Z[i, j + 2 * n] <- 0
      Z[j + 2 * n, i] <- 0
      
      Z[j, i + 2 * n] <- 0 
      Z[i + 2 * n, j] <- 0 
      
      Z[i + n, j + 2 * n] <- 0
      Z[j + 2 * n, i + n] <- 0
      
      Z[j + n, i + 2 * n] <- 0
      Z[i + 2 * n, j + n] <- 0
    }
  }
}

for (i in 1:n){
  Z[i, i + n] <- exp(-T_scale)
  Z[i + n, i] <- exp(-T_scale)
  
  Z[i + 2 * n, i] <- exp(-2 * T_scale)
  Z[i, i + 2 * n] <- exp(-2 * T_scale)
  
  Z[i + n, i + 2 * n] <- exp(-T_scale)
  Z[i + 2 * n, i + n] <- exp(-T_scale)
}
for (i in 1:(3 * n)){
  S[i] <- sum(Z[i, ])
}
for (i in 1:(3 * n)){
  for (j in 1:(3 * n)){
    W_star[i, j] <- Z[i, j]/S[i]
  }
}


      
      # read data: train with NA
      train_DF_full <- subset(data_train_DF, year == 2015 & month == 9)
      # covariate scaling
      covariate <- names(train_DF_full)
      for (i in 8:length(covariate)){
        ori <- train_DF_full[, i]
        scaled <- (ori - mean(ori))/sd(ori)
        train_DF_full[, i] <- scaled
      }
      
      
      train_DF_pre_origin <- subset(data_train_DF, year == 2011 & month == 9)
      train_DF_pre <- subset(train_DF_pre_origin, select = -c(CNT, lon, lat, year, month))
      train_DF_later_origin <- subset(data_train_DF, year == 2013 & month == 9)
      train_DF_later <- subset(train_DF_later_origin, select = -c(CNT, lon, lat, year, month))
      
      
      # merge year.temp - 2, year.temp, year.temp - 1 (merge three years)
      train_DF <- subset(train_DF_full, select = -c(CNT, lon, lat, year, month))
      train_DF_merge <- rbind(train_DF_pre, train_DF, train_DF_later)
          
  
  
  
    ## BA prediction
    Y <- train_DF_merge$BA
    Z <- log(1 + Y)
    train_DF_origin <- train_DF_new <- train_DF_merge
    train_DF_merge <- rename(train_DF_merge, "WZ" = "BA")
    train_DF_merge$WZ <- W_star %*% Z
    tau <- 0.9
    X <- as.matrix(train_DF_merge[, -1])
    instru <- as.data.frame(W_star %*% X)  # instrument variable
    names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                       "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                       "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                       "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
    train_with_instru <- cbind(train_DF_merge, instru)
    fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru)
    WZ_hat <- fit_WZ$fitted.values
    if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 => Z(tau) = 0
      pred[, j] <- rep(0, nrow(test_DF))
      next
    }
    train_DF_new <- rename(train_DF_new, "Z" = "BA")
    train_DF_new$Z <- Z
    train_DF_new$WZ_hat <- WZ_hat
    fit_Z <- rq(Z ~., tau = tau, data = train_DF_new)
    summary(fit_Z)
    
  

    
    
    

    ################# CNT
    
    # read data: train with NA
    train_DF_full <- subset(data_train_DF, year == 2015 & month == 9)
    # covariate scaling
    covariate <- names(train_DF_full)
    for (i in 8:length(covariate)){
      ori <- train_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF_full[, i] <- scaled
    }
    
    
    train_DF_pre_origin <- subset(data_train_DF, year == 2011 & month == 9)
    train_DF_pre <- subset(train_DF_pre_origin, select = -c(BA, lon, lat, year, month))
    train_DF_later_origin <- subset(data_train_DF, year == 2013 & month == 9)
    train_DF_later <- subset(train_DF_later_origin, select = -c(BA, lon, lat, year, month))
    
    
    # merge 2011, 2015, 2013(merge three years)
    train_DF <- subset(train_DF_full, select = -c(BA, lon, lat, year, month))
    train_DF_merge <- rbind(train_DF_pre, train_DF, train_DF_later)
    
    
    
    
    ## CNT prediction
    Y <- train_DF_merge$CNT
    Z <- Y + runif(length(Y), 0, 0.999)
    train_DF_origin <- train_DF_new <- train_DF_merge
    train_DF_merge <- rename(train_DF_merge, "WZ" = "CNT")
    train_DF_merge$WZ <- W_star %*% Z
    tau <- 0.5
    X <- as.matrix(train_DF_merge[, -1])
    instru <- as.data.frame(W_star %*% X)  # instrument variable
    names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                       "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                       "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                       "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
    train_with_instru <- cbind(train_DF_merge, instru)
    fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "sfn")
    WZ_hat <- fit_WZ$fitted.values
    if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 => Z(tau) = 0
      pred[, j] <- rep(0, nrow(test_DF))
      next
    }
    train_DF_new <- rename(train_DF_new, "Z" = "CNT")
    train_DF_new$Z <- Z
    train_DF_new$WZ_hat <- WZ_hat
    fit_Z <- rq(Z ~., tau = tau, data = train_DF_new)
    summary(fit_Z)
    
    
    
    