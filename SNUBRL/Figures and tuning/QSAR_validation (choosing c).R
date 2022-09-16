## Validation CNT choosing c.
#rm(list=ls())

library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
#library(spdep)
#library(caret)
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




setwd("~/Desktop")
load("data_full.RData")  # full data without missing

#load("data_train.RData")
#is.na((data_train_DF))

data_ori <- data_DF
train_temp <- subset(data_ori, year == 1993 & month == 3)
loc <- cbind(train_temp$lon, train_temp$lat)
n <- nrow(loc)

# construct weight matrix (cut-off value is 800km)
W <- z <- matrix(rep(0, n * n), nr = n, nc = n)
s <- c()
sigma <- 270     # scale 270 (tuning parameter) 
for (i in 2:n){  
  for (j in 1:(i - 1)){
    dist <- distGeo(loc[i, ], loc[j, ])/1000  # distance between ith and jth locations. (km)
    if (dist < 800){
      z[i, j] <- exp(-dist/sigma)   
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






## Validation of CNT

#setwd("~/Desktop")
load("data_train.RData")   # data_train_DF (full data with NA's)
data_odd <- readRDS("data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year %% 2 == 1)


#### Tuning parameter for extreme value index gamma
c <- 0.1 * (1:50)
#c.temp <- c[15] 
#print(c.temp)
####

validation_c <- c()


# validation c
for (index in c(11, 12, 13, 14, 16, 17, 18, 19,21,22,23,24,26,27,28,29,31,32,33,34,36,37,38,39,41,42,43,44,46,47,48,49)){
c.temp <- c[index]

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
    train <- subset(data_odd, year == year.temp & month == month.temp)
    train <- subset(train, select = -c(lon, lat, year, month))
    
    # covariate scaling
    covariate <- names(train)
    
    for (i in 3:length(covariate)){
      ori <- train[, i]
      train[, i] <- (ori - mean(ori))/sd(ori)
    }
    train_ori <- train
    missing <- is.na(train$CNT)
    
    # BA=0 index
    BA_0 <- (train[missing, ]$BA == 0)  
    BA_0[is.na(BA_0)] <- 0
    
    loc_obs <- loc[!missing, ]
    loc_miss <- loc[missing, ]
    y <- train_ori$CNT[!missing]
    y <- as.matrix(y)
    train <- train_ori[!missing, ]
    
    # true data 
    test <- subset(data_odd_true, year == year.temp & month == month.temp)
    test <- subset(test, select = -c(BA, lon, lat, year, month))
    true_cnt <- test$CNT
    # covariate scaling
    covariate <- names(test)
    
    for (i in 2:length(covariate)){
      ori <- test[, i]
      test[, i] <- (ori - mean(ori))/sd(ori)
    }
    test <- test_ori <- test[missing, ]
    
    
    
    ## missing location (Thin-plate spline, Tps)
    y <- train$CNT
    y <- as.matrix(y)
    #Poissonization (Anscombe transform)
    z <- sqrt(y + 3/8)
    fit <- Tps(loc_obs, z, lambda = 0, lon.lat = TRUE, miles = FALSE)
    temp <- predict(fit, loc_miss)
    #impute <- floor(temp)
    impute <- floor(temp^2 - 1/8)
    #summary(impute)
    #surface(fit)
    #US(add = TRUE, col = "grey")
    for (i in 1:nrow(impute)){ # If impute value is negative, set 0.
      if (impute[i] < 0){
        impute[i] <- 0
      }
    }
    train_imp <- train_ori
    train_imp$CNT[missing] <- impute
    
    # train with imputation
    train <- train_imp
    train <- subset(train, select = -BA)
    # test set without CNT
    test <- subset(test, select = -CNT)
    
    ## CNT prediction
    Y <- train$CNT
    
    # jittering (discrete -> continuous)
    Z <- Y + runif(length(Y), 0, 0.999)
    # Poissonization (Anscombe transform)
    #Z <- sqrt(Y + 3/8)
    #summary(Z)
    taus <- seq(0.05, 0.95, by = 0.01)
    m <- length(taus)
    pred <- matrix(nr = nrow(test), nc = m)   # estimated quantiles
    train_ori <- train_new <- train
    train <- rename(train, "WZ" = "CNT")
    train$WY <- W %*% Z
    for (j in 1:m){
      #j <- 1
      #j <- j + 1
      tau <- taus[j]
      #print(taus[j])
      train <- train_new <- train_ori
      train <- rename(train, "WZ" = "CNT")
      train$WZ <- W %*% Z
      X <- as.matrix(train[, -1])
      instru <- as.data.frame(W %*% X)  # instrument variable
      names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                         "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                         "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                         "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
      train_with_instru <- cbind(train, instru)
      fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "sfn")
      WZ_hat <- fit_WZ$fitted.values
      if (sum(WZ_hat^2) < 1e-10){   # WY_hat(tau) = 0 implies Y(tau) = 0
        pred[, j] <- rep(0, nrow(test))
        next
      }
      train_new <- rename(train_new, "Z" = "CNT")
      train_new$Z <- Z
      train_new$WZ_hat <- WZ_hat
      fit_Z <- rq(Z ~., tau = tau, data = train_new, method = "sfn")
      test$WZ_hat <- WZ_hat[missing]
      result <- predict(fit_Z, test)
      pred[, j] <- result
    }
    
   
    
    for (i in 1:nrow(test)){  # sort quantiles
      pred[i, ] <- sort(pred[i, ])
    }
    
    for (i in 1:nrow(test)){  # negative quantile set to be 0.
      for (j in 1:length(taus)){
        if (pred[i, j] < 0){
          pred[i, j] <- 0
        }
      }
    }
    
    
    # Extreme value index (gamma) estimation 
    gamma <- c()
    for (i in 1:3503){ # ith location
      #i <- 945
      obs.temp <- subset(data_odd, lon == loc[i, 1] & lat == loc[i, 2])
      sort.temp <- sort(obs.temp$CNT, decreasing = TRUE)
      T.temp <- length(sort.temp)
      k.temp <- floor(c.temp * T.temp^{1/3})
      gamma.temp <- 0
      if (sort.temp[k.temp] < 1e-10){
        gamma[i] <- 0 
      } else {
        for (t in 1:k.temp){
          gamma.temp <- gamma.temp + log(sort.temp[t]/sort.temp[k.temp])
        }
        gamma[i] <- gamma.temp/k.temp
      }
    }
    # gamma
    gamma_missing <- gamma[missing]
    
    
    
    # Extremal conditional quantile extrapolation
    # gamma <- 0.1   # extremal index (to be estimated by location)
    
    tau_0 <- 0.95
    taus_ext <- c(0.96, 0.97, 0.98, 0.99)
    pred.ori <- pred
    for (j in 1:length(taus_ext)){
      tau_ext_temp <- taus_ext[j]
      multipli <- ((1 - tau_0)/(1 - tau_ext_temp))^(gamma_missing) 
      pred <- cbind(pred, pred.ori[, length(taus)] * multipli) # extrapolate 0.96, 0.97, 0.98, 0.99th quantiles
    }
    
    pred.ori <- pred
    for (j in 1:length(taus_ext)){
      tau_ext_temp <- taus_ext[j]
      multipli <- ((1 - tau_0)/(1 - tau_ext_temp))^(gamma_missing) 
      pred <- cbind(pred.ori[, 1] / multipli, pred) # extrapolate 0.96, 0.97, 0.98, 0.99th quantiles
    }
    
    taus <- c(sort(1 - taus_ext), taus, taus_ext)
    tau_0 <- 0.99
    taus_ext2 <- c(0.995, 0.999, 0.9995)
    
    pred.ori <- pred
    for (j in 1:length(taus_ext2)){
      tau_ext_temp <- taus_ext2[j]
      multipli <- ((1 - tau_0)/(1 - tau_ext_temp))^(gamma_missing) 
      pred <- cbind(pred, pred.ori[, length(taus)] * multipli) # extrapolate 0.995, 0.999, and 0.9995th quantiles
    }
    
    pred.ori <- pred
    for (j in 1:length(taus_ext2)){
      tau_ext_temp <- taus_ext2[j]
      multipli <- ((1 - tau_0)/(1 - tau_ext_temp))^(gamma_missing) 
      pred <- cbind(pred.ori[, 1]/multipli, pred) # extrapolate 0.995, 0.999, and 0.9995th quantiles
    }
    taus <- c(0, sort(1 - taus_ext2), taus, taus_ext2)
    pred <- cbind(rep(0, nrow(test)), pred)   # add tau = 0
    
    
    for (i in 1:nrow(pred)){ # find quantile of CNT (=Y)
      pred[i, ] <- floor(pred[i, ])
    }
    
    for (i in 1:nrow(test)){  # negative value is set to be 0.
      for (j in 1:length(taus)){
        if (pred[i, j] < 0){
          pred[i, j] <- 0
        }
      }
    }
    
    
    
    
    # calculate score (CNT)
    score_cnt <- 0
    for (i in 1:nrow(test)){
      #print(i)
      for (k in 1:length(u_cnt)){
        #print(k)
        if (as.logical(BA_0)[i]){ 
          p_hat <- 1 
        } else if (pred[i, length(taus)] < u_cnt[k]){
          p_hat <- 1
        } else {
          p_hat <- taus[p_hat_find(pred[i, ], u_cnt[k])]
        }
        #print(score_cnt)
        score_cnt <- score_cnt + weights_cnt[k] * (as.numeric(true_cnt[i] <= u_cnt[k]) - p_hat)^2 
      }
    }
    #score_cnt
  
    cnt.score[iter] <- score_cnt
    print(sum(cnt.score[1:iter]))
    
    iter <- iter + 1
  }
}

validation_c[index] <- sum(cnt.score)

}

print(validation_c)

