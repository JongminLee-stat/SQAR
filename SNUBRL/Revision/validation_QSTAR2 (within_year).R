## validation (QSTAR2, within a year)

rm(list=ls())
library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(patchwork)
library(splines2)
library(fields)
#library(caret)


p_hat_find <- function (vec, b){
  vec_rev <- vector(length = length(vec))
  for (i in 1:length(vec)){
    vec_rev[i] <- vec[length(vec) + 1 - i]
  }
  ind <- length(vec) - which.min(abs(vec_rev - b)) + 1
  return(ind)
}

# modify your working directory with the data file if necessary (using the setwd() command):

#setwd("/Users/jongminlee/Desktop/EVA2021/Data")
load("data_train.RData")
data_odd_DF <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)
train_DF_full <- subset(data_odd_DF, year == 1995 & month == 5)
loc <- cbind(train_DF_full$lon, train_DF_full$lat)
n <- nrow(loc)


# weight matrix (cut-off value is 800km)
# d <- matrix(rep(0, n * n), nr = n, nc = n)
W_star <- Z <- matrix(rep(0, (6 * n)^2), nr = 6 * n, nc = 6 * n)
S <- c()
scale <- 270    # (distance) scale hyperparameter
T <- 30         # time scale hyperparameter  
year_T <- 7 * T
T_scale <- T/scale  # scaled time distance 
year_T_scale <- year_T/scale
for (i in 2:n){  
  for (j in 1:(i - 1)){
    dist <- distGeo(loc[i, ], loc[j, ])/1000  # distance between i-th and j-th locations (km)
    dist_scale <- dist/scale
    if (dist < sqrt(4800^2 + (year_T_scale + 2 * T_scale)^2)/6){   # About 1/6 of the largest distance set to be cutoff value. 
      rep0 <- exp(-dist_scale)
      Z[i, j] <- rep0
      Z[j, i] <- rep0
      Z[i + n, j + n] <- rep0
      Z[j + n, i + n] <- rep0
      Z[i + 2 * n, j + 2 * n] <- rep0
      Z[j + 2 * n, i + 2 * n] <- rep0
      Z[i + 3 * n, j + 3 * n] <- rep0
      Z[j + 3 * n, i + 3 * n] <- rep0
      Z[i + 4 * n, j + 4 * n] <- rep0
      Z[j + 4 * n, i + 4 * n] <- rep0
      Z[i + 5 * n, j + 5 * n] <- rep0
      Z[j + 5 * n, i + 5 * n] <- rep0
      
      rep1 <- exp(-sqrt(T_scale^2 + dist_scale^2)) 
      Z[i, j + n] <- rep1
      Z[j + n, i] <- rep1 
      Z[j, i + n] <- rep1   
      Z[i + n, j] <- rep1   
      
      rep2 <- exp(-sqrt(4 * T_scale^2 + dist_scale^2))
      Z[i, j + 2 * n] <- rep2 
      Z[j + 2 * n, i] <- rep2
      Z[j, i + 2 * n] <- rep2 
      Z[i + 2 * n, j] <- rep2 
      
      Z[i + n, j + 2 * n] <- rep1
      Z[j + 2 * n, i + n] <- rep1
      Z[j + n, i + 2 * n] <- rep1
      Z[i + 2 * n, j + n] <- rep1
      
      rep3 <- exp(-sqrt(year_T_scale^2 + dist_scale^2))
      Z[j, i + 3 * n] <- rep3 
      Z[i + 3 * n, j] <- rep3
      Z[i, j + 3 * n] <- rep3 
      Z[j + 3 * n, i] <- rep3
      
      rep4 <- exp(-sqrt((year_T_scale + T_scale)^2 + dist_scale^2))
      Z[j + 3 * n, i + n] <- rep4  
      Z[i + n, j + 3 * n] <- rep4
      Z[i + 3 * n, j + n] <- rep4
      Z[j + n, i + 3 * n] <- rep4
      
      rep5 <- exp(-sqrt((year_T_scale + 2 * T_scale)^2 + dist_scale^2))
      Z[j + 2 * n, i + 3 * n] <- rep5
      Z[i + 3 * n, j + 2 * n] <- rep5
      Z[i + 2 * n, j + 3 * n] <- rep5
      Z[j + 3 * n, i + 2 * n] <- rep5
      
      Z[j, i + 4 * n] <- rep4  
      Z[i + 4 * n, j] <- rep4
      Z[i, j + 4 * n] <- rep4 
      Z[j + 4 * n, i] <- rep4
      
      Z[j + 4 * n, i + n] <- rep3 
      Z[i + n, j + 4 * n] <- rep3
      Z[i + 4 * n, j + n] <- rep3 
      Z[j + n, i + 4 * n] <- rep3
      
      Z[j + 4 * n, i + 2 * n] <- rep4 
      Z[i + 2 * n, j + 4 * n] <- rep4
      Z[i + 4 * n, j + 2 * n] <- rep4 
      Z[j + 2 * n, i + 4 * n] <- rep4
      
      Z[j + 4 * n, i + 3 * n] <- rep1
      Z[i + 3 * n, j + 4 * n] <- rep1
      Z[i + 4 * n, j + 3 * n] <- rep1
      Z[j + 3 * n, i + 4 * n] <- rep1
      
      Z[j, i + 5 * n] <- rep5 
      Z[i + 5 * n, j] <- rep5
      Z[i, j + 5 * n] <- rep5 
      Z[j + 5 * n, i] <- rep5
      
      Z[j + 5 * n, i + n] <- rep4 
      Z[i + n, j + 5 * n] <- rep4
      Z[i + 5 * n, j + n] <- rep4 
      Z[j + n, i + 5 * n] <- rep4
      
      Z[j + 5 * n, i + 2 * n] <- rep3 
      Z[i + 2 * n, j + 5 * n] <- rep3
      Z[i + 5 * n, j + 2 * n] <- rep3 
      Z[j + 2 * n, i + 5 * n] <- rep3
      
      Z[j + 5 * n, i + 3 * n] <- rep2
      Z[i + 3 * n, j + 5 * n] <- rep2
      Z[i + 5 * n, j + 3 * n] <- rep2
      Z[j + 3 * n, i + 5 * n] <- rep2
      
      Z[j + 4 * n, i + 5 * n] <- rep1
      Z[i + 5 * n, j + 4 * n] <- rep1
      Z[i + 4 * n, j + 5 * n] <- rep1
      Z[j + 5 * n, i + 4 * n] <- rep1
      
    }
  }
}

for (i in 1:n){
  rep_1 <- exp(-T_scale)
  Z[i, i + n] <- rep_1
  Z[i + n, i] <- rep_1
  Z[i + n, i + 2 * n] <- rep_1
  Z[i + 2 * n, i + n] <- rep_1
  Z[i + 3 * n, i + 4 * n] <- rep_1
  Z[i + 4 * n, i + 3 * n] <- rep_1
  Z[i + 4 * n, i + 5 * n] <- rep_1
  Z[i + 5 * n, i + 4 * n] <- rep_1

  rep_2 <- rep_1^2
  Z[i, i + 2 * n]  <- rep_2
  Z[i + 2 * n, i]  <- rep_2
  Z[i + 3 * n, i + 5 * n]  <- rep_2
  Z[i + 5 * n, i + 3 * n]  <- rep_2
  
  rep_3 <- exp(-sqrt((year_T_scale + T_scale)^2 + dist_scale^2))
  Z[i + n, i + 3 * n]  <- rep_3
  Z[i + 3 * n, i + n]  <- rep_3
  Z[i + 4 * n, i + 2 * n]  <- rep_3
  Z[i + 2 * n, i + 4 * n]  <- rep_3
  
  rep_4 <- exp(-sqrt((year_T_scale + 2 * T_scale)^2 + dist_scale^2))
  Z[i + 2 * n, i + 3 * n] <- rep_4
  Z[i + 3 * n, i + 2 * n] <- rep_4
  
  rep_5 <- exp(-sqrt(year_T_scale^2 + dist_scale^2))
  Z[i + 3 * n, i] <- rep_5
  Z[i, i + 3 * n] <- rep_5
  Z[i + 4 * n, i + n] <- rep_5
  Z[i + n, i + 4 * n] <- rep_5
  Z[i + 5 * n, i + 2 * n] <- rep_5
  Z[i + 2 * n, i + 5 * n] <- rep_5

  Z[i + 4 * n, i] <- rep_3
  Z[i, i + 4 * n] <- rep_3
  Z[i + 5 * n, i + n] <- rep_3
  Z[i + n, i + 5 * n] <- rep_3

  Z[i + 5 * n, i] <- rep_4 
  Z[i, i + 5 * n] <- rep_4 
}
for (i in 1:(6 * n)){
  S[i] <- sum(Z[i, ])
}
for (i in 1:(6 * n)){
  for (j in 1:(6 * n)){
    W_star[i, j] <- Z[i, j]/S[i]
  }
}



year <- c(1993, 1997, 2001, 2005, 2009, 2013)  # estimated years
month <- seq(3, 9, by = 1)
t <- length(year) * length(month) 
ba.score <- vector(length = t)
iter <- 1
for (year_temp in year){
  for (month_temp in month){
    print(iter)
    #year_temp <- 1993
    #month_temp <- 3
    if (year_temp == 1993){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else { # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 1997){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2001){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2005){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2009){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else {  # year_temp = 2013
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else { # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    }
    
    
    # covariate scaling
    covariate <- names(train_DF)
    for (i in 8:length(covariate)){
      ori <- train_DF[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF[, i] <- scaled
    }
    missing <- is.na(train_DF$BA)
    
    # CNT=0 index
    CNT_0 <- (train_DF[missing, ]$CNT == 0)  
    CNT_0[is.na(CNT_0)] <- 0
    
    # imputation
    impute <- c()
    for (j in 1:3){
      missing_temp <- missing[(3503 * (j - 1) + 1):(3503 * j)]
      loc_obs <- loc[!missing_temp, ]
      loc_miss <- loc[missing_temp, ]
      # imputation at missing location by month (Thin-plate spline, Tps)
      y <- train_DF[c(rep(FALSE, 3503 * (j - 1)), !missing_temp, rep(FALSE, 3503 * (6 - j))), ]$BA
      y <- as.matrix(y)
      y_log <- log(1 + y)
      fit <- Tps(loc_obs, y_log, lambda = 0, lon.lat = TRUE, miles = FALSE)
      impute_log <- predict(fit, loc_miss)
      impute_temp <- exp(impute_log) - 1
      for (i in 1:nrow(impute_temp)){  # If imputed value is negative, set 0.
        if (impute_temp[i] < 0){
           impute_temp[i] <- 0
        }
      }
      impute <- c(impute, as.vector(impute_temp))
    }
    #length(impute)
    
    train_DF$BA[missing] <- impute 
    train_DF <- subset(train_DF, select = -c(CNT, lon, lat, year, month))
    
    if (month_temp == 3){
      missing_ori <- c(missing[1:3503], rep(FALSE, 3503 * 5))
    } else if (month_temp == 9){
      missing_ori <- c(rep(FALSE, 3503 * 2), missing[(3503 * 2 + 1):(3503 * 3)], rep(FALSE, 3503 * 3))
    } else {
      missing_ori <- c(rep(FALSE, 3503), missing[(3503 + 1):(3503 * 2)], rep(FALSE, 3503 * 4))
    }
    
    # true data 
    test_DF <- subset(data_odd_true, year == year_temp & month == month_temp)
    
    # covariate scaling
    covariate <- names(test_DF)
    for (i in 8:length(covariate)){
      ori <- test_DF[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      test_DF[, i] <- scaled
    }
    
    test_DF <- subset(test_DF, select = -c(CNT, lon, lat, year, month))
    if (month_temp == 3){
      test_DF <- test_DF_origin <- test_DF[missing[1:3503], ]
    } else if (month_temp == 9){
      test_DF <- test_DF_origin <- test_DF[missing[(3503 * 2 + 1):(3503 * 3)], ]
    } else {
      test_DF <- test_DF_origin <- test_DF[missing[(3503 + 1):(3503 * 2)], ]
    }
    # test set without BA
    test_DF <- subset(test_DF, select = -BA)
    
    
    ## BA prediction
    Y <- train_DF$BA
    Z <- log(1 + Y)
    taus <- c(seq(0.05, 0.95, by = 0.05))
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "BA")
    train_DF$WZ <- W_star %*% Z
    
    for (j in 1:m){
      tau <- taus[j]
      print(taus[j])
      train_DF <- train_DF_new <- train_DF_origin
      train_DF <- rename(train_DF, "WZ" = "BA")
      train_DF$WZ <- W_star %*% Z
      X <- as.matrix(train_DF[, -1])
      instru <- as.data.frame(W_star %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train_DF, instru)
      fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "sfn")
      WZ_hat <- fit_WZ$fitted.values
      if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 => Z(tau) = 0
        pred[, j] <- rep(0, nrow(test_DF))
        next
      }
      train_DF_new <- rename(train_DF_new, "Z" = "BA")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      
      fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
      test_DF$WZ_hat <- WZ_hat[missing_ori]
      
      result <- predict(fit_Z, test_DF)
      pred[, j] <- result
    }
    for (i in 1:nrow(test_DF)){  # sort quantiles
      pred[i, ] <- sort(pred[i, ])
    }
    for (i in 1:nrow(test_DF)){  # negative quantile is set to be 0.
      for (j in 1:length(taus)){
        if (pred[i, j] < 0){
          pred[i, j] <- 0
        }
      }
    }
    
    
    # extremal quantile extrapolation
    gamma <- 0.1   # extremal index (to be estimated)
    tau_0 <- 0.95
    taus_ext <- c(0.96, 0.97, 0.98, 0.99)
    multipli <- ((1 - tau_0)/(1 - taus_ext))^(gamma)
    pred.old <- pred
    for (i in 1:length(taus_ext)){  # extrapolate 0.96, 0.97, 0.98, 0.99
      pred <- cbind(pred, pred.old[, length(taus)] * multipli[i])
    }
    taus <- c(taus, taus_ext)
    tau_0 <- 0.99
    taus_ext2 <- c(0.995, 0.999, 0.9995)
    multipli2 <- ((1 - tau_0)/(1 - taus_ext2))^(gamma) 
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.995, 0.999, 0.9995
      pred <- cbind(pred, pred.old[, length(taus)] * multipli2[i])
    }
    taus <- c(0, taus, taus_ext2)
    pred <- cbind(rep(0, nrow(test_DF)), pred)   # add tau=0
  
    # calculate score (BA)
    score_ba <- 0
    for (i in 1:nrow(test_DF)){
      for (k in 1:length(u_ba)){
        if (as.logical(CNT_0)[i]){ 
          p_hat <- 1 
        } else if ( pred[i, length(taus)] < log(1 + u_ba[k]) ){
          p_hat <- 1
        } else {
          p_hat <- taus[p_hat_find(pred[i, ], log(1 + u_ba[k]))]
        }
        score_ba <- score_ba + weights_ba[k] * (as.numeric(test_DF_origin$BA[i] <= u_ba[k]) - p_hat)^2 
      }
    }
    
    ba.score[iter] <- score_ba
    iter <- iter + 1
  }
 }

##
print("ba.score is")  # scale=270, T=30
sum(ba.score)









#### CNT prediction

# modify your working directory with the data file if necessary (using the setwd() command):
#setwd("/Users/jongminlee/Desktop/EVA2021/Data")
load("data_train.RData")
data_odd_DF <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)
train_DF_full <- subset(data_odd_DF, year == 1993 & month == 3)
loc <- cbind(train_DF_full$lon, train_DF_full$lat)
n <- nrow(loc)

year <- c(1993, 1997, 2001, 2005, 2009, 2013)
month <- seq(3, 9, by = 1)
t <- length(year)
cnt.score <- vector(length = t)
iter <- 1
for (year_temp in year){
  for (month_temp in month){
    #year_temp <- 1993
    #month_temp <- 3
    print(iter)
    if (year_temp == 1993){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else { # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 1997){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2001){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 2) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2005){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp + 10) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else if (year_temp == 2009){
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else {   # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 6) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    } else {  # year_temp = 2013
      if (month_temp == 3){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= (month_temp + 2))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & month <= (month_temp + 2))
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else if (month_temp == 9){
        train_DF_1 <- filter(data_odd_DF, year == year_temp & month <= month_temp & (month_temp - 2) <= month)
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & month <= month_temp & (month_temp - 2) <= month)
        train_DF <- rbind(train_DF_1, train_DF_2)
      } else { # month_temp = 4, 5, 6, 7, 8
        train_DF_1 <- filter(data_odd_DF, year == year_temp & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF_2 <- filter(data_odd_DF, year == (year_temp - 10) & (month_temp - 1) <= month & month <= (month_temp + 1))
        train_DF <- rbind(train_DF_1, train_DF_2)
      }
    }
    
    
    # covariate scaling
    covariate <- names(train_DF)
    for (i in 8:length(covariate)){
      ori <- train_DF[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF[, i] <- scaled
    }
    missing <- is.na(train_DF$CNT)
    
    # BA=0 index
    BA_0 <- (train_DF[missing, ]$BA == 0)  
    BA_0[is.na(BA_0)] <- 0
    
    # imputation
    impute <- c()
    for (j in 1:3){
      missing_temp <- missing[(3503 * (j - 1) + 1):(3503 * j)]
      loc_obs <- loc[!missing_temp, ]
      loc_miss <- loc[missing_temp, ]
      # imputation at missing location by month (Thin-plate spline, Tps)
      y <- train_DF[c(rep(FALSE, 3503 * (j - 1)), !missing_temp, rep(FALSE, 3503 * (6 - j))), ]$CNT
      y <- as.matrix(y)
      z <- sqrt(y + 3/8)  # Poissonization (Anscombe transform)
      fit <- Tps(loc_obs, z, lambda = 0, lon.lat = TRUE, miles = FALSE)
      temp <- predict(fit, loc_miss)
      impute_temp <- floor(temp^2 - 1/8)
      for (i in 1:nrow(impute_temp)){  # If impute value is negative, set 0.
        if (impute_temp[i] < 0){
          impute_temp[i] <- 0
        }
      }
      impute <- c(impute, as.vector(impute_temp))
    }
    # train with imputation
    train_DF$CNT[missing] <- impute
    train_DF <- subset(train_DF, select = -c(BA, lon, lat, year, month))
    if (month_temp == 3){
      missing_ori <- c(missing[1:3503], rep(FALSE, 3503 * 5))
    } else if (month_temp == 9){
      missing_ori <- c(rep(FALSE, 3503 * 2), missing[(3503 * 2 + 1):(3503 * 3)], rep(FALSE, 3503 * 3))
    } else {
      missing_ori <- c(rep(FALSE, 3503), missing[(3503 + 1):(3503 * 2)], rep(FALSE, 3503 * 4))
    }
    
    # true data 
    test_DF <- subset(data_odd_true, year == year_temp & month == month_temp)

    # covariate scaling
    covariate <- names(test_DF)
    for (i in 8:length(covariate)){
      ori <- test_DF[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      test_DF[, i] <- scaled
    }
    test_DF <- subset(test_DF, select = -c(BA, lon, lat, year, month))
    if (month_temp == 3){
      test_DF <- test_DF_origin <- test_DF[missing[1:3503], ]
    } else if (month_temp == 9){
      test_DF <- test_DF_origin <- test_DF[missing[(3503 * 2 + 1):(3503 * 3)], ]
    } else {
      test_DF <- test_DF_origin <- test_DF[missing[(3503 + 1):(3503 * 2)], ]
    }
    
    # test set without CNT
    test_DF <- subset(test_DF, select = -CNT)
    
    ## CNT prediction
    Y <- train_DF$CNT
    # Poissonization (Anscombe transform)
    Z <- sqrt(Y + 3/8)
    taus <- seq(0.05, 0.95, by = 0.05)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "CNT")
    train_DF$WZ <- W_star %*% Z
    for (j in 1:m){
      tau <- taus[j]
      print(taus[j])
      train_DF <- train_DF_new <- train_DF_origin
      train_DF <- rename(train_DF, "WZ" = "CNT")
      train_DF$WZ <- W_star %*% Z
      X <- as.matrix(train_DF[, -1])
      instru <- as.data.frame(W_star %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train_DF, instru)
      fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "sfn")
      
      WZ_hat <- fit_WZ$fitted.values
      if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 implies Z(tau) = 0
        pred[, j] <- rep(0, nrow(test_DF))
        next
      }
      train_DF_new <- rename(train_DF_new, "Z" = "CNT")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
      test_DF$WZ_hat <- WZ_hat[missing_ori]
      result <- predict(fit_Z, test_DF)
      pred[, j] <- result
    }
    
    for (i in 1:nrow(test_DF)){  # sort quantiles
      pred[i, ] <- sort(pred[i, ])
    }
    
    for (i in 1:nrow(test_DF)){  # negative quantile is set to be 0.
      for (j in 1:length(taus)){
        if (pred[i, j] < 0){
          pred[i, j] <- 0
        }
      }
    }
    
    
    # extremal quantile extrapolation
    gamma <- 0.1   # extremal index (to be estimated)
    tau_0 <- 0.95
    taus_ext <- c(0.96, 0.97, 0.98, 0.99)
    multipli <- ((1 - tau_0)/(1 - taus_ext))^(gamma)
    pred.old <- pred
    for (i in 1:length(taus_ext)){  # extrapolate 0.96, 0.97, 0.98, 0.99
      pred <- cbind(pred, pred.old[, length(taus)] * multipli[i])
    }
    taus <- c(taus, taus_ext)
    tau_0 <- 0.99
    taus_ext2 <- c(0.995, 0.999, 0.9995)
    multipli2 <- ((1 - tau_0)/(1 - taus_ext2))^(gamma) 
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.995, 0.999, 0.9995
      pred <- cbind(pred, pred.old[, length(taus)] * multipli2[i])
    }
    taus <- c(0, taus, taus_ext2)
    pred <- cbind(rep(0, nrow(test_DF)), pred)   # add tau=0
    
    for (i in 1:nrow(pred)){ # find quantile of CNT (=Y)
      pred[i, ] <- floor((pred[i, ])^2 - 1/8)    # use -1/8, rather than -3/8
    }
    
    
    # calculate score (CNT)
    score_cnt <- 0
    for (i in 1:nrow(test_DF)){
      for (k in 1:length(u_cnt)){
        if (as.logical(BA_0)[i]){ 
          p_hat <- 1 
        } else if (pred[i, length(taus)] < u_cnt[k]){
          p_hat <- 1
        } else {
          p_hat <- taus[p_hat_find(pred[i, ], u_cnt[k])]
        }
        score_cnt <- score_cnt + weights_cnt[k] * (as.numeric(test_DF_origin$CNT[i] <= u_cnt[k]) - p_hat)^2 
      }
    }
    cnt.score[iter] <- score_cnt
    iter <- iter + 1
  }  
}

print("cnt score is")
sum(cnt.score)   # scale=270, T=30


