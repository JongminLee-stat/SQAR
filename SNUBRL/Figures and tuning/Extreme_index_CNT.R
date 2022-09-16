## Extreme value index (gamma) estimate (CNT)

rm(list=ls())

library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
#library(caret)
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
train <- subset(data_ori, year == 1993 & month == 3)
loc <- cbind(train$lon, train$lat)
train <- subset(data_ori, year == year.temp & month == month.temp)
train <- subset(train, select = -c(BA, lon, lat, year, month))



c.temp <- c <- 0.6
# Extreme value index (gamma) estimation 
gamma <- c()
for (i in 1:3503){ # ith location
  #i <- 945
  obs.temp <- subset(data_ori, lon == loc[i, 1] & lat == loc[i, 2])
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

    

    
    
 #pdf("gamma_CNT.pdf", height = 20, width = 15)
 map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
 map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
 quilt.plot(loc, gamma, ny = 35, add = TRUE)
 title("Extreme value index (CNT)", cex.main = 2)
     
     
     
    
   