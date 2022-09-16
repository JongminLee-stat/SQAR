## validation (QSTAR model, seasnality)

rm(list=ls())


library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(caret)
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
scale <- 100    # (distance) scale 
T <- 0            
T_scale <- T/scale  # scaled time distance, (2 year) 
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
#for (i in 2:(3 * n)){
#  for (j in 1:(i - 1)){
#    Z[j, i] <- Z[i, j]
#  }
#}

for (i in 1:(3 * n)){
  S[i] <- sum(Z[i, ])
}

for (i in 1:(3 * n)){
  for (j in 1:(3 * n)){
    W_star[i, j] <- Z[i, j]/S[i]
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
    
   
    # missing location (Thin-plate spline, Tps)
    y <- train_DF$BA
    y <- as.matrix(y)
    y_log <- log(1 + y)
    fit <- Tps(loc_obs, y_log, lambda = 0, lon.lat = TRUE, miles = FALSE)
    impute_log <- predict(fit, loc_miss)
    impute <- exp(impute_log) - 1
    #summary(impute)
    #surface(fit)
    #US(add = TRUE, col = "grey")
    for (i in 1:nrow(impute)){ # If impute value is negative, set 0.
      if (impute[i] < 0){
        impute[i] <- 0
      }
    }
    train_imp <- train_DF_full
    train_imp$BA[missing] <- impute
    
    # train with imputation
    train_DF <- train_imp
    train_DF <- subset(train_DF, select = -c(CNT, lon, lat, year, month))
    # test set without BA
    test_DF <- subset(test_DF, select = -BA)
    
    if (year.temp != 1993){
      train_DF_pre_origin <- subset(data_odd_DF, year == 2015 & month == month.temp)
      train_DF_pre <- subset(train_DF_pre_origin, select = -c(CNT, lon, lat, year, month))
    } else {
      train_DF_pre_origin <- subset(data_odd_DF, year == year.temp + 2 & month == month.temp)
      train_DF_pre <- subset(train_DF_pre_origin, select = -c(CNT, lon, lat, year, month))
    }
    
    train_DF_later_origin <- subset(data_odd_DF, year == year.temp + 2 & month == month.temp)
    train_DF_later <- subset(train_DF_later_origin, select = -c(CNT, lon, lat, year, month))
    
    # merge year.temp - 2, year.temp, year.temp + 2 (merge three years)
    train_DF_merge <- rbind(train_DF_pre, train_DF, train_DF_later)
    
    
    ## BA prediction
    Y <- train_DF_merge$BA
    Z <- log(1 + Y)
    
    
    taus <- seq(0.05, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF_merge
    train_DF_merge <- rename(train_DF_merge, "WZ" = "BA")
    train_DF_merge$WZ <- W_star %*% Z
    
    
    for (j in 1:m){
      #j <- 1
      tau <- taus[j]
      #print(taus[j])
      train_DF_merge <- train_DF_new <- train_DF_origin
      train_DF_merge <- rename(train_DF_merge, "WZ" = "BA")
      head(train_DF_merge)
      train_DF_merge$WZ <- W_star %*% Z
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
      train_DF_new <- rename(train_DF_new, "Z" = "BA")
      train_DF_new$Z <- Z
      train_DF_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
      missing_merge <- c(rep(FALSE, 3503), missing, rep(FALSE, 3503))
      test_DF$WZ_hat <- WZ_hat[missing_merge]
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
    pred.old <- pred
    for (i in 1:length(taus_ext)){  # extrapolate 0.04, 0.03, 0.02, 0.01
      pred <- cbind(pred.old[, 1]/multipli[i], pred)
    }
    
    taus <- c(sort(1 - taus_ext), taus, taus_ext)
    tau_0 <- 0.99
    taus_ext2 <- c(0.995, 0.999, 0.9995)
    multipli2 <- ((1 - tau_0)/(1 - taus_ext2))^(gamma) 
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.995, 0.999, 0.9995
      pred <- cbind(pred, pred.old[, length(taus)] * multipli2[i])
    }
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.001, 0.005, 0.0005
      pred <- cbind(pred.old[, 1]/multipli2[i], pred)
    }
    
    
    taus <- c(0, sort(1 - taus_ext2), taus, taus_ext2)
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
print("ba.score is")  # scale = 100
sum(ba.score)







#### CNT prediction

# modify your working directory with the data file if necessary (using the setwd() command):

#getwd()
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
    
    
    ## missing location (Thin-plate spline, an imputation method for count is needed. shirink transformation of CNT is also needed.)
    y <- train_DF$CNT
    y <- as.matrix(y)
    # Poissonization (Anscombe transform)
    z <-  sqrt(y + 3/8)
    fit <- Tps(loc_obs, z, lambda = 0, lon.lat = TRUE, miles = FALSE)
    temp <- predict(fit, loc_miss)
    impute <- floor(temp^2 - 1/8)  # use -1/8 rather than -3/8 (unbiased)
    #summary(impute)
    #surface(fit)
    #US(add = TRUE, col = "grey")
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
    
    if (year.temp != 1993){
      train_DF_pre_origin <- subset(data_odd_DF, year == 2015 & month == month.temp)
      train_DF_pre <- subset(train_DF_pre_origin, select = -c(BA, lon, lat, year, month))
    } else {
      train_DF_pre_origin <- subset(data_odd_DF, year == year.temp + 2 & month == month.temp)
      train_DF_pre <- subset(train_DF_pre_origin, select = -c(BA, lon, lat, year, month))
    }
    
    train_DF_later_origin <- subset(data_odd_DF, year == year.temp + 2 & month == month.temp)
    train_DF_later <- subset(train_DF_later_origin, select = -c(BA, lon, lat, year, month))
    
    # merge year.temp - 2, year.temp, year.temp + 2 (merge three years)
    train_DF_merge <- rbind(train_DF_pre, train_DF, train_DF_later)
    
    
    
    
    
    ## CNT prediction
    Y <- train_DF_merge$CNT
    Z <- sqrt(Y + 3/8)                           # Poissonization (Anscombe transform)      
    taus <- seq(0.05, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test_DF), nc = m)   # estimated quantiles
    train_DF_origin <- train_DF_new <- train_DF_merge
    train_DF_merge <- rename(train_DF_merge, "WZ" = "CNT")
    train_DF_merge$WZ <- W_star %*% Z
    for (j in 1:m){
      tau <- taus[j]
      #print(taus[j])
      train_DF_merge <- train_DF_new <- train_DF_origin
      train_DF_merge <- rename(train_DF_merge, "WZ" = "CNT")
      train_DF_merge$WZ <- W_star %*% Z
      X <- as.matrix(train_DF_merge[, -1])
      instru <- as.data.frame(W_star %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train_DF_merge, instru)
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
      missing_merge <- c(rep(FALSE, 3503), missing, rep(FALSE, 3503))
      test_DF$WZ_hat <- WZ_hat[missing_merge]
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
    pred.old <- pred
    for (i in 1:length(taus_ext)){  # extrapolate 0.04, 0.03, 0.02, 0.01
      pred <- cbind(pred.old[, 1]/multipli[i], pred)
    }
    
    taus <- c(sort(1 - taus_ext), taus, taus_ext)
    tau_0 <- 0.99
    taus_ext2 <- c(0.995, 0.999, 0.9995)
    multipli2 <- ((1 - tau_0)/(1 - taus_ext2))^(gamma) 
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.995, 0.999, 0.9995
      pred <- cbind(pred, pred.old[, length(taus)] * multipli2[i])
    }
    pred.old <- pred
    for (i in 1:length(taus_ext2)){ # extrapolate 0.001, 0.005, 0.0005
      pred <- cbind(pred.old[, 1]/multipli2[i], pred)
    }
    
    taus <- c(0, sort(1 - taus_ext2), taus, taus_ext2)
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
sum(cnt.score)   # scale=100


