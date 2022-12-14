## validation (QSAR estimated by IVQR method)

rm(list=ls())
library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(patchwork)
library(splines2)
library(spam)
library(fields)
library(magrittr)




p_hat_find <- function (vec, b){
  vec_rev <- vector(length = length(vec))
  for (i in 1:length(vec)){
    vec_rev[i] <- vec[length(vec) + 1 - i]
  }
  ind <- length(vec) - which.min(abs(vec_rev - b)) + 1
  return(ind)
}

# modify your working directory with the data file if necessary (using the setwd() command):
getwd()
#setwd("C:/Users/KIAS/Desktop")
load("data_train.RData")
data_odd_DF <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)
train_DF_full <- subset(data_odd_DF, year == 1993 & month == 3)
loc <- cbind(train_DF_full$lon, train_DF_full$lat)
n <- nrow(loc)



# construct weight matrix (cut-off value is 800km)
W <- z <- matrix(rep(0, n * n), nr = n, nc = n)
s <- c()
scale <- 270     # scale 270 (tuning hyperparameter) 
for (i in 2:n){  
  for (j in 1:(i - 1)){
    dist <- distGeo(loc[i, ], loc[j, ])/1000  # distance between ith and jth locations. (km)
    if (dist < 800){
      z[i, j] <- exp(-dist/scale)   
    } else {
      z[i, j] <- 0
    }
    z[j, i] <- z[i, j]
  }
}
for (i in 1:n){
  s[i] <- sum(z[i, ])
}
for (i in 1:n){
  for (j in 1:n){
    W[i, j] <- z[i, j]/s[i]
  }
}
year <- c(1993, 1997, 2001, 2005, 2009, 2013)
month <- seq(3, 9, by = 1)
t <- length(year) * length(month)
ba.score <- vector(length = t)
iter <- 1

for (year.temp in year){
  for (month.temp in month){
    #year.temp <- 1993
    #month.temp <- 3
    print(iter)
    
    # read data: train with NA
    train_DF_full <- subset(data_odd_DF, year == year.temp & month == month.temp)
    
    # covariate scaling
    covariate <- names(train_DF_full)
    for (i in 8:length(covariate)){
      ori <- train_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF_full[, i] <- scaled
    }
    missing <- is.na(train_DF_full$BA)
 
    # CNT = 0 index
    CNT_0 <- (train_DF_full[missing, ]$CNT == 0)  
    CNT_0[is.na(CNT_0)] <- 0
    loc_obs <- loc[!missing, ]
    loc_miss <- loc[missing, ]
    y <- train_DF_full$BA[!missing]
    y <- as.matrix(y)
    train_DF <- train_DF_full[!missing, ]
  
    # true data 
    test_DF_full <- subset(data_odd_true, year == year.temp & month == month.temp)
    # covariate scaling
    covariate <- names(test_DF_full)
    for (i in 8:length(covariate)){
      ori <- test_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      test_DF_full[, i] <- scaled
    }
    test_DF_full <- subset(test_DF_full, select = -c(CNT, lon, lat, year, month))
    test_DF <- test_DF_origin <- test_DF_full[missing, ]
    
    # missing location (Thin-plate spline, Tps)
    y <- train_DF$BA
    y <- as.matrix(y)
    y_log <- log(1 + y)
    fit <- Tps(loc_obs, y_log, lambda = 0, lon.lat = TRUE, miles = FALSE)
    impute_log <- predict(fit, loc_miss)
    impute <- exp(impute_log) - 1
    for (i in 1:nrow(impute)){ # If impute value is negative, set 0.
      if (impute[i] < 0){
        impute[i] <- 0
      }
    }
    train_imp <- train_DF_full
    train_imp$BA[missing] <- impute
    train_DF <- train_imp                     # train with imputation
    train_DF <- subset(train_DF, select = -c(CNT, lon, lat, year, month))
    test_DF <- subset(test_DF, select = -BA)  # test set without BA
    
    ## BA prediction
    Y <- train_DF$BA
    Z <- log(1 + Y)
    taus <- seq(0.01, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "BA")
    train_DF$WZ <- W %*% Z
    
    
    for (j in 1:m){
      #j <- 1
      tau <- taus[j]
      print(taus[j])
      lambda <- seq(0.10, 0.90, by = 0.05)  # grid search (candidates of lambda)
      memory_temp <- vector(length=length(lambda))
      iter2 <- 0
      for (lambda_temp in lambda){ # Do QR (Z - lambdaWZ ~ X, WX)
        iter2 <- iter2 + 1  
        #lambda_temp <- 0.5  
        #print(lambda_temp)
        train_DF <- train_DF_origin
        train_DF <- rename(train_DF, "first" = "BA")
        train_DF$first <- (Z - lambda_temp * W %*% Z)
        X <- as.matrix(train_DF[, -1])
        instru <- as.data.frame(W %*% X)  # instrumental variable
        names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
        train_with_instru <- cbind(train_DF, instru)
        fit_first <- rq(first ~., tau = tau, data = train_with_instru, method = "sfn")
        memory_temp[iter2] <- sqrt(sum((fit_first$coefficients[33:63])^2))
      }
        lambda_hat <- lambda[which.min(memory_temp)]   # IVQR estimate of lambda
        train_with_instru <- rename(train_with_instru, "second" = "first")
        train_with_instru$second <- (Z - lambda_hat * W %*% Z)
        fit_IVQR <- rq(second ~., tau = tau, data = train_with_instru, method = "sfn")
        #beta_hat <- fit_IVQR$coefficients[2:32]  # estimate of beta
        
        # estimation of WY_hat
        train_DF <- train_DF_origin
        train_DF <- rename(train_DF, "WZ" = "BA")
        train_DF$WZ <- W %*% Z
        X <- as.matrix(train_DF[, -1])
        instru <- as.data.frame(W %*% X)  # instrument variable
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
        pred[, j] <- (lambda_hat * WZ_hat[missing] + predict(fit_IVQR, train_with_instru)[missing]) 
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
    
    # Extrapolate the extremal quantile based on intermediate quantile
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
    pred <- cbind(rep(0, nrow(test_DF)), pred)   # add tau = 0
    
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
    print("year_temp")
    print("month_temp")
    print("score_ba is")
    print(score_ba)
  }
 }

print("Validation score of BA (IVQR) is")
sum(ba.score)





## CNT prediction (IVQR estimation) 

# modify your working directory with the data file if necessary (using the setwd() command):
load("data_train.RData")
data_odd_DF <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)
train_DF_full <- subset(data_odd_DF, year == 1993 & month == 3)
loc <- cbind(train_DF_full$lon, train_DF_full$lat)
n <- nrow(loc)


year <- c(1993, 1997, 2001, 2005, 2009, 2013)
month <- seq(3, 9, by = 1)
t <- length(year) * length(month)
cnt.score <- vector(length = t)
iter <- 1

for (year.temp in year){
  for (month.temp in month){
    #print(year.temp)
    #print(month.temp)
    print(iter)
    
    # read data: train with NA
    train_DF_full <- subset(data_odd_DF, year == year.temp & month == month.temp)
    # covariate scaling
    covariate <- names(train_DF_full)
    for (i in 8:length(covariate)){
      ori <- train_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF_full[, i] <- scaled
    }
    missing <- is.na(train_DF_full$CNT)
    
    # BA=0 index
    BA_0 <- (train_DF_full[missing, ]$BA == 0)  
    BA_0[is.na(BA_0)] <- 0
    loc_obs <- loc[!missing, ]
    loc_miss <- loc[missing, ]
    y <- train_DF_full$CNT[!missing]
    y <- as.matrix(y)
    train_DF <- train_DF_full[!missing, ]
    
    # true data 
    test_DF_full <- subset(data_odd_true, year == year.temp & month == month.temp)
    # covariate scaling
    covariate <- names(test_DF_full)
    for (i in 8:length(covariate)){
      ori <- test_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      test_DF_full[, i] <- scaled
    }
    test_DF_full <- subset(test_DF_full, select = -c(BA, lon, lat, year, month))
    test_DF <- test_DF_origin <- test_DF_full[missing, ]
    
    
    ## missing location (Tps)
    y <- train_DF$CNT
    y <- as.matrix(y)
    # Poissonization (Anscombe transform)
    z <- sqrt(y + 3/8)
    fit <- Tps(loc_obs, z, lambda = 0, lon.lat = TRUE, miles = FALSE)
    temp <- predict(fit, loc_miss)
    impute <- floor(temp^2 - 1/8)
    for (i in 1:nrow(impute)){ # If impute value is negative, set 0.
      if (impute[i] < 0){
        impute[i] <- 0
      }
    }
    train_imp <- train_DF_full
    train_imp$CNT[missing] <- impute
    
    # train with imputation
    train_DF <- train_imp
    train_DF <- subset(train_DF, select = -c(BA, lon, lat, year, month))
    # test set without CNT
    test_DF <- subset(test_DF, select = -CNT)
    
    ## CNT prediction
    Y <- train_DF$CNT
    # jittering (discrete -> continuous)
    Z <- Y + runif(length(Y), 0, 0.999)
    #Z <- sqrt(Y + 3/8) # Poissonization (Anscombe transform)
    taus <- seq(0.01, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "CNT")
    train_DF$WY <- W %*% Z
    for (j in 1:m){
      tau <- taus[j]
      #print(iter)
      print(taus[j])
      lambda <- seq(0.10, 1.10, by = 0.05)  # grid search (candidates of lambda)
      memory_temp <- vector(length=length(lambda))
      iter2 <- 0
      for (lambda_temp in lambda){ # Do quantile regression (Z - lambdaWZ ~ X, WX)
        iter2 <- iter2 + 1  
        train_DF <- train_DF_origin
        train_DF <- rename(train_DF, "first" = "CNT")
        train_DF$first <- (Z - lambda_temp * W %*% Z)
        X <- as.matrix(train_DF[, -1])
        instru <- as.data.frame(W %*% X)  # instrumental variable
        names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                           "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                           "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                           "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
        train_with_instru <- cbind(train_DF, instru)
        fit_first <- rq(first ~., tau = tau, data = train_with_instru, method = "sfn")
        memory_temp[iter2] <- sqrt(sum((fit_first$coefficients[33:63])^2))
      }
      lambda_hat <- lambda[which.min(memory_temp)]   # IVQR estimate of lambda
      train_with_instru <- rename(train_with_instru, "second" = "first")
      train_with_instru$second <- (Z - lambda_hat * W %*% Z)
      fit_IVQR <- rq(second ~., tau = tau, data = train_with_instru, method = "sfn")
      
      # estimation of WY_hat
      train_DF <- train_DF_origin
      train_DF <- rename(train_DF, "WZ" = "CNT")
      train_DF$WZ <- W %*% Z
      X <- as.matrix(train_DF[, -1])
      instru <- as.data.frame(W %*% X)  # instrument variable
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
      pred[, j] <- (lambda_hat * WZ_hat[missing] + predict(fit_IVQR, train_with_instru)[missing]) 
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

    # extremal conditional quantile extrapolation
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
    pred <- cbind(rep(0, nrow(test_DF)), pred)   # add tau = 0
 
    for (i in 1:nrow(pred)){ # find quantile of CNT (=Y)
      pred[i, ] <- floor(pred[i, ])
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
    print("year.temp")
    print("month.temp")
    print("score_cnt is")
    print(score_cnt)
    cnt.score[iter] <- score_cnt
    iter <- iter + 1
  }
}


#### Score of CNT
print("Validation score (IVQR) of CNT is")
sum(cnt.score)


#### Score of BA
print("Validation score (IVQR) of BA is")
sum(ba.score)


