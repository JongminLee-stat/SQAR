## Estimation of 99.95 quantile response CNT

#rm(list=ls())
library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
#library(spdep)
library(caret)
library(patchwork)
library(splines2)
library(fields)
library(maps)

n <- 3503
year.temp <- 2015
month.temp <- 9
# modify your working directory with the data file if necessary (using the setwd() command):

load("data_full.RData")  # full data without missing
data_ori <- data_DF
train <- subset(data_ori, year == year.temp & month == month.temp)
loc <- cbind(train$lon, train$lat)
train <- subset(data_ori, year == year.temp & month == month.temp)
train <- subset(train, select = -c(year, month, BA, lon, lat))

## construct weight matrix (cut-off value is 800km)
W <- z <- matrix(rep(0, n * n), nr = n, nc = n)
s <- c()
scale <- 270     # scale 270 (tuning parameter) 
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



tau <- 0.95 
gamma <- 0.1 
Z <- (train$CNT + runif(length(train$CNT), 0, 0.999))

train_DF <- train_DF_new <- train
train_DF <- rename(train_DF, "WZ" = "CNT")
train_DF$WZ <- W %*% Z
X <- as.matrix(train_DF[, -1])
instru <- as.data.frame(W %*% X)  # instrumental variable
names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                   "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                   "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                   "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
train_with_instru <- cbind(train_DF, instru)
fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "sfn")
WZ_hat <- fit_WZ$fitted.values
train_DF_new <- rename(train_DF_new, "Z" = "CNT")
train_DF_new$Z <- Z
train_DF_new$WZ_hat <- WZ_hat
fit_Z <- rq(Z ~., tau = tau, data = train_DF_new, method = "sfn")
train <- subset(train, select = -CNT)
train$WZ_hat <- WZ_hat
Z <- predict(fit_Z, train)
Z_CNT <- floor(Z * ((1-0.95)/(1-0.99))^(gamma) * ((1-0.99)/(1-0.9995))^(gamma))



# BA
#par(mfrow = c(1, 1),
#    mar = c(5.5, 1.5, 5.5, 1.5)
#) 
map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
quilt.plot(loc, Z_CNT, ny = 35, add = TRUE)
title("99.95% quantiles of jittered CNT", cex.main = 1.6)


