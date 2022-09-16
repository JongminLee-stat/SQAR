## validation (Iterative SQAR model)

rm(list=ls())

library(imputeTS)
library(SparseM)
library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(patchwork)
library(splines2)
library(fields)

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
    print(iter)
    #year.temp <- 1993
    #month.temp <- 3
    
    train_DF_full <- subset(data_odd_DF, year == year.temp & month == month.temp)
    
    # covariate scaling
    covariate <- names(train_DF_full)
    for (i in 8:length(covariate)){
      ori <- train_DF_full[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_DF_full[, i] <- scaled
    }
    missing <- is.na(train_DF_full$BA)
    
    # CNT=0 index
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
    
    #ts_temp <- c()
    #impute <- vector(length = nrow(loc_miss))
    #for (i in nrow(loc_miss)){
    #  ts_temp <- subset(data_odd_DF, lon == loc_miss[i, 1] & lat == loc_miss[i, 2])$BA
    #  log_ts_temp <- log(1 + ts_temp)
    #  impute[i] <- exp(na_interpolation(log_ts_temp, option = "spline")[((year.temp - 1993)/2) * 7 + month.temp - 2]) - 1
    #  if (impute[i] < 0){
    #    impute[i] <- 0
    #  }
    #}
  
    # missing location (Thin-plate spline, Tps)
    y <- train_DF$BA
    y <- as.matrix(y)
    y_log <- log(1 + y)
    impute_new <- impute_log <- impute <- vector(length = nrow(loc_miss))
    fit <- Tps(loc_obs, y_log, lambda = 0, lon.lat = TRUE, miles = FALSE)
    impute_log <- predict(fit, loc_miss)
    impute <- exp(impute_log) - 1
    for (i in 1:nrow(impute)){ # If imputed value is negative, set 0.
      if (impute[i] < 0){
        impute[i] <- 0
      }
    }
    # train with imputation
    train_DF <- subset(train_DF_full, select = -c(CNT, lon, lat, year, month))
    train_DF$BA[missing] <- as.vector(impute)
    # test set without BA
    test_DF <- subset(test_DF, select = -BA)
    
    
    
    ## BA prediction
    Y <- train_DF$BA
    Z <- log(1 + Y)
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "BA")
    train_DF$WZ <- W %*% Z

    ## repeative impute
    for (i in 1:3){
      train_DF <- train_DF_new <- train_DF_origin
      # head(train_DF)
      Y <- train_DF$BA
      Z <- log(1 + Y)
      train_DF <- rename(train_DF, "WZ" = "BA")
      train_DF$WZ <- W %*% Z
      X <- as.matrix(train_DF[, -1])
      instru <- as.data.frame(W %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train_DF, instru)
      fit_WZ <- rq(WZ ~., tau = 0.5, data = train_with_instru, method = "sfn")
      WZ_hat <- fit_WZ$fitted.values
      if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 => Z(tau) = 0
        pred[, j] <- rep(0, nrow(test_DF))
        next
      }
      train_DF_new <- rename(train_DF_new, "Z" = "BA")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = 0.5, data = train_DF_new, method = "sfn")
      test_DF$WZ_hat <- WZ_hat[missing]
      impute_log <- predict(fit_Z, test_DF)
      impute <- exp(impute_log) - 1
      for (i in 1:length(impute)){ # If imputed value is negative, set 0.
        if (impute[i] < 0){
          impute[i] <- 0
        }
      }
      
      train_DF$WZ[missing] <- impute
      train_DF <- rename(train_DF, "BA" = "WZ")
      #train_DF <- subset(train_DF, select = -c(CNT, lon, lat, year, month))
      train_DF_origin <- train_DF_new <- train_DF
    }
    
    
    
    ## BA prediction
    Y <- train_DF$BA
    Z <- log(1 + Y)
    taus <- seq(0.01, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    
    
    for (j in 1:m){
      tau <- taus[j]
      #print(taus[j])
      train_DF <- train_DF_new <- train_DF_origin
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
      train_DF_new <- rename(train_DF_new, "Z" = "BA")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
      
      test_DF$WZ_hat <- WZ_hat[missing]
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
    score_ba
    ba.score[iter] <- score_ba
    iter <- iter + 1
  }
}

print("Validation score of BA is")
sum(ba.score)




## CNT prediction

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
t <- length(year) * length(month)
cnt.score <- vector(length = t)
iter <- 1

for (year.temp in year){
  for (month.temp in month){
    print(iter)
    #year.temp <- 1993
    #month.temp <- 3
    
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
    
    #ts_temp <- c()
    #impute <- vector(length = nrow(loc_miss))
    #for (i in nrow(loc_miss)){
    #  ts_temp <- subset(data_odd_DF, lon == loc_miss[i, 1] & lat == loc_miss[i, 2])$BA
    #  log_ts_temp <- log(1 + ts_temp)
    #  impute[i] <- exp(na_interpolation(log_ts_temp, option = "spline")[((year.temp - 1993)/2) * 7 + month.temp - 2]) - 1
    #  if (impute[i] < 0){
    #    impute[i] <- 0
    #  }
    #}
    
    ## missing location (Thin-plate spline, Tps)
    y <- train_DF$CNT
    y <- as.matrix(y)
    z <- sqrt(y + 3/8) # Poissonization (Anscombe transform)
    fit <- Tps(loc_obs, z, lambda = 0, lon.lat = TRUE, miles = FALSE)
    temp <- predict(fit, loc_miss)
    impute <- floor(temp^2 - 1/8)
    
    for (i in 1:nrow(impute)){ # If imputed value is negative, set 0.
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
    #Z <- Y + runif(length(Y), 0, 0.999)
    Z <- sqrt(Y + 3/8) # Poissonization (Anscombe transform)
    taus <- seq(0.01, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF
    train_DF <- rename(train_DF, "WZ" = "CNT")
    train_DF$WY <- W %*% Z
    
    
    ## repeative impute
    for (i in 1:3){
      train_DF <- train_DF_new <- train_DF_origin
      # head(train_DF)
      Y <- train_DF$CNT
      Z <- sqrt(Y + 3/8)
      train_DF <- rename(train_DF, "WZ" = "CNT")
      train_DF$WZ <- W %*% Z
      X <- as.matrix(train_DF[, -1])
      instru <- as.data.frame(W %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train_DF, instru)
      fit_WZ <- rq(WZ ~., tau = 0.5, data = train_with_instru, method = "sfn")
      WZ_hat <- fit_WZ$fitted.values
      if (sum(WZ_hat^2) < 1e-10){   # WZ_hat(tau) = 0 => Z(tau) = 0
        pred[, j] <- rep(0, nrow(test_DF))
        next
      }
      train_DF_new <- rename(train_DF_new, "Z" = "CNT")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = 0.5, data = train_DF_new, method = "sfn")
      test_DF$WZ_hat <- WZ_hat[missing]
      impute_log <- predict(fit_Z, test_DF)
      impute <- exp(impute_log) - 1
      for (i in 1:length(impute)){ # If imputed value is negative, set 0.
        if (impute[i] < 0){
          impute[i] <- 0
        }
      }
      
      train_DF$WZ[missing] <- impute
      train_DF <- rename(train_DF, "CNT" = "WZ")
      #train_DF <- subset(train_DF, select = -c(BA, lon, lat, year, month))
      train_DF_origin <- train_DF_new <- train_DF
    }
    
    
    ## CNT prediction
    Y <- train_DF$CNT
    Z <- sqrt(Y + 3/8)
    taus <- seq(0.01, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    for (j in 1:m){
      tau <- taus[j]
      #print(taus[j])
      train_DF <- train_DF_new <- train_DF_origin
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
      if (sum(WZ_hat^2) < 1e-10){   # WY_hat(tau) = 0 implies Y(tau) = 0
        pred[, j] <- rep(0, nrow(test_DF))
        next
      }
      train_DF_new <- rename(train_DF_new, "Z" = "CNT")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
      test_DF$WZ_hat <- WZ_hat[missing]
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
    pred <- cbind(rep(0, nrow(test_DF)), pred)   # add tau=0
 
    for (i in 1:nrow(pred)){ # find quantile of CNT (=Y)
      pred[i, ] <- floor((pred[i, ])^2 - 1/8)
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

print("Validation score of CNT is")
sum(cnt.score)


