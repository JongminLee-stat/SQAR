## EDA (quilt plot) (Quantile SAR model)
rm(list=ls())

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
#setwd()


#### Uncertainty plot of BA (western)
load("data_full.RData")  # full data without missing
year <- seq(1993, 2015, by = 1)
month <- seq(3, 9, by = 1)
Z_BA <- vector(length = 161)
iter <- 0
for (year.temp in year){
  for(month.temp in month){  
    #year.temp <- 1993
    #month.temp <- 3
    iter <- iter + 1
    data_ori <- data_DF
    train <- subset(data_ori, year == year.temp & month == month.temp)
    select <- subset(train, lon == -123.25 & lat == 38.75)
    Z_BA[iter] <- log(1 + select$BA)
  }
}


boxplot(Z_BA, cex.lab = 1.3, cex.main = 1.25, ylab = "log(1+BA)", main = "Boxplot of log-transformed BA (western)")



#### Uncertainty plot of BA (central)
load("data_full.RData")  # full data without missing
year <- seq(1993, 2015, by = 1)
month <- seq(3, 9, by = 1)
Z_BA <- vector(length = 161)
iter <- 0
for (year.temp in year){
  for(month.temp in month){  
    #year.temp <- 1993
    #month.temp <- 3
    iter <- iter + 1
    data_ori <- data_DF
    train <- subset(data_ori, year == year.temp & month == month.temp)
    select <- subset(train, lon == -95.75 & lat == 38.75)
    Z_BA[iter] <- log(1 + select$BA)
    
  }
}

boxplot(Z_BA, cex.lab = 1.3, cex.main = 1.25, ylab = "log(1+BA)", main = "Boxplot of log-transformed BA (central)")

