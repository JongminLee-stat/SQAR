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

# Western region
Z_BA1 <- vector(length = 161)
iter <- 0
for (year.temp in year){
  for(month.temp in month){  
    #year.temp <- 1993
    #month.temp <- 3
    iter <- iter + 1
    data_ori <- data_DF
    train <- subset(data_ori, year == year.temp & month == month.temp)
    select <- subset(train, lon == -123.25 & lat == 38.75)
    Z_BA1[iter] <- log(1 + select$BA)
  }
}

# Central region
year <- seq(1993, 2015, by = 1)
month <- seq(3, 9, by = 1)
Z_BA2 <- vector(length = 161)
iter <- 0
for (year.temp in year){
  for(month.temp in month){  
    #year.temp <- 1993
    #month.temp <- 3
    iter <- iter + 1
    data_ori <- data_DF
    train <- subset(data_ori, year == year.temp & month == month.temp)
    select <- subset(train, lon == -95.75 & lat == 38.75)
    Z_BA2[iter] <- log(1 + select$BA)
    
  }
}

BA_df <- data.frame(Western = Z_BA1, Central = Z_BA2)

pdf("uncertainty_BA.pdf", height=5, width=5)
boxplot(BA_df, cex.lab = 1.6, cex.main = 1.5, cex.axis = 1.4,
        ylab = "log(1+BA)", main = "Boxplot of log-transformed BA")
dev.off()
