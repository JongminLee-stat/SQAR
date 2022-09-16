## Boxplots of Spatial lag parameter (BA) (tau=0.5, 0.9)
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

W <- z <- matrix(rep(0, n * n), nr = n, nc = n)
s <- c()
sigma <- 270     # sigma 270 (tuning parameter) 
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



year <- 1993:2015
month <- 3:9
n <- 3503
# modify your working directory with the data file if necessary (using the setwd() command):
load("data_full.RData")  # full data without missing

data_ori <- data_DF
lambda <- c()
iter <- 1 
for(year.temp in year){
  for(month.temp in month){
    #year.temp <- 2015    
    #month.temp <- 9
    print(iter)
    train <- subset(data_ori, year == year.temp & month == month.temp)
    train <- subset(train, select = -c(CNT, lon, lat, year, month))
    train_new <- train
    
    ## prediction of distribution of BA 
    Y <- train$BA
    Z <- log(1 + Y)
    tau <- 0.5
    
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
    lambda[iter] <- as.numeric((fit_Z$coefficients)[33])
    
    iter <- iter + 1 
  }
}

lambda3 <- lambda[seq(1, 155, by = 7)]  # March
lambda4 <- lambda[seq(2, 156, by = 7)]  # April
lambda5 <- lambda[seq(3, 157, by = 7)]  # May
lambda6 <-lambda[seq(4, 158, by = 7)]   # June
lambda7 <- lambda[seq(5, 159, by = 7)]  # July
lambda8 <-lambda[seq(6, 160, by = 7)]   # August
lambda9 <-lambda[seq(7, 161, by = 7)]   # September

## plot by month
plot(x = 1993:2015, y = lambda3, xlab = "March", main = "Estimates of spatial lag paramter (March)", type = "o")
plot(x = 1993:2015, y = lambda4, xlab = "April", main = "Estimates of spatial lag paramter (April)", type = "o")
plot(x = 1993:2015, y = lambda5, xlab = "May", main = "Estimates of spatial lag paramter (May)", type = "o")
plot(x = 1993:2015, y = lambda6, xlab = "June", main = "Estimates of spatial lag paramter (June)", type = "o")
plot(x = 1993:2015, y = lambda7, xlab = "July", main = "Estimates of spatial lag paramter (July)", type = "o")
plot(x = 1993:2015, y = lambda8, xlab = "August", main = "Estimates of spatial lag paramter (August)", type ="o")
plot(x = 1993:2015, y = lambda9, xlab = "September", main = "Estimates of spatial lag paramter (September)", type = "o")



#### Boxplot by month
pdf("lambda_month_BA.pdf", height = 10, width = 15)
par(mfrow = c(1, 1),
    mar = c(5.5, 5.5, 5.5, 5.5)
)        
boxplot(lambda3, lambda4, lambda5, lambda6, lambda7, lambda8, lambda9, 
        xlab = "Month", axes = FALSE, ylab = "Lambda", cex.lab = 3, cex.main = 3, main = "Spatial lag parameter (tau=0.5) by month (BA)")
axis(1, at = 1:7, cex.axis = 2.5, lab = c("3", "4", "5", "6", "7", "8", '9'))
axis(2, at = 0.05 * 0:20, cex.axis = 2.5)
dev.off()


#### Boxplot by month
pdf("lambda_year_BA.pdf", height = 10, width = 15)
par(mfrow = c(1, 1),
    mar = c(5.5, 5.5, 5.5, 5.5)
)        
boxplot(lambda[1:7], lambda[8:14],lambda[15:21], lambda[22:28], lambda[29:35],lambda[36:42],lambda[43:49],lambda[50:56],lambda[57:63],lambda[64:70],
        lambda[71:77],lambda[78:84],lambda[85:91],lambda[92:98],lambda[99:105],lambda[106:112],lambda[113:119],lambda[120:126],
        lambda[127:133], lambda[134:140], lambda[141:147],lambda[148:154],lambda[155:161], names = 1993:2015, xlab = "Year", ylab = "Lambda",
        axes = FALSE, cex.lab = 3, cex.main = 3, main = "Spatial lag parameter (tau=0.5) by year (BA)")
axis(1, at = 1:23, cex.axis = 2.5, lab = c("93", "94", "95", "96", "97", "98", "99", "00", "01", "02", "03", 
                                           "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15"))
axis(2, at = 0.05 * 0:20, cex.axis = 2.5)
dev.off()

