## Coefficient varying (BA)
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




# modify your working directory with the data file if necessary (using the setwd() command):
load("data_full.RData")  # full data without missing
data_ori <- data_DF
train <- subset(data_ori, year == 1993 & month == 3)

loc <- cbind(train$lon, train$lat)
n <- 3503
# Construct weight matrix (cut-off value is 800km)
#W <- z <- matrix(rep(0, n * n), nr = n, nc = n)
#s <- c()
#sigma <- 270     # sigma 270 (tuning parameter) 
#for (i in 2:n){  
#  for (j in 1:(i - 1)){
#    dist <- distGeo(loc[i, ], loc[j, ])/1000  # distance between ith and jth locations. (km)
#    if (dist < 800){
#      z[i, j] <- exp(-dist/sigma)   
#    } else {
#      z[i, j] <- 0
#    }
#    z[j, i] <- z[i, j]
#  }
#}
#for (i in 1:n){
#  s[i] <- sum(z[i, ])
#}
#for (i in 1:n){
#  for (j in 1:n){
#    W[i, j] <- z[i, j]/s[i]
#  }
#}


year.temp <- 2015    
month.temp <- 9
data_ori <- data_DF
coeff <- c()
iter <- 1 
taus <- 1:99 * 0.01

for(tau.temp in taus){
  
    train <- subset(data_ori, year == year.temp & month == month.temp)
    train <- subset(train, select = -c(CNT, lon, lat, year, month))
    ## prediction of distribution of BA 
    Y <- train$BA
    Z <- log(1 + Y)
  
    tau <- tau.temp  # quantile level of interest
    #print(iter)
    
    train_ori <- train_new <- train
    train <- rename(train, "WZ" = "BA")
    # covariate scaling
    covariate <- names(train)
    for (i in 1:length(covariate)){
      ori <- train[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train[, i] <- scaled 
    }
    
    train$WZ <- W %*% Z
    X <- as.matrix(train[, -1])
    instru <- as.data.frame(W %*% X)  # instrument variable
    names(instru) <- c("area.inst", "lc1.inst", "lc2.inst", "lc3.inst", "lc4.inst", "lc5.inst", "lc6.inst", "lc7.inst", "lc8.inst", "lc9.inst",
                       "lc10.inst", "lc11.inst", "lc12.inst", "lc13.inst", "lc14.inst", "lc15.inst", "lc16.inst", "lc17.inst", "lc18.inst",
                       "altiMean.inst", "altiSD.inst", "clim1.inst", "clim2.inst", "clim3.inst", "clim4.inst", "clim5.inst", "clim6.inst",
                       "clim7.inst", "clim8.inst", "clim9.inst", "clim10.inst")
    train_with_instru <- cbind(train, instru)
    fit_WZ <- rq(WZ ~., tau = tau, data = train_with_instru, method = "fn")
    WZ_hat <- fit_WZ$fitted.values
    train_new <- rename(train, "Z" = "WZ")
    train_new$Z <- Z
    train_new$WZ_hat <- WZ_hat
    
    # covariate scaling
    covariate <- names(train_new)
    for (i in 1:length(covariate)){
      ori <- train_new[, i]
      scaled <- (ori - mean(ori))/sd(ori)
      train_new[, i] <- scaled 
    }
    
    fit_Z <- rq(Z ~., tau = tau, data = train_new, method = "fn")
    fit_Z$coefficients[33]
    
    coeff[iter] <- as.numeric((fit_Z$coefficients)[33])
    iter <- iter + 1 
  }

pdf("coeff_lambda_BA.pdf", height = 10, width = 15)
par(mfrow = c(1, 1),
    mar = c(6, 6, 6, 6)
)   
plot(x = 38:99 * 0.01, y = coeff[38:99], xlim = c(0.38, 1), xlab = "Quantile level", ylab = "Coefficient", main = "Spatial lagged dependent variable (BA)", 
     cex.axis = 2, cex.lab = 2.4, cex.main = 3.5, type = "l", lty = "solid")
dev.off()


