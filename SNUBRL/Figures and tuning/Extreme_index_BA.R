## Extreme value index (gamma) estimate (CNT)

rm(list=ls())

#install.packages("geosphere")
#install.packages("quantreg")
#install.packages("dplyr")
#install.packages("splines2")
#install.packages("patchwork")
#install.packages("spdep")
#install.packages("purrr")


library(geosphere)
library(quantreg)
library(dplyr)
library(psych)
library(spdep)
library(caret)
library(patchwork)
library(splines2)
library(fields)
library(maps)
library(purrr)

n <- 3503
year.temp <- 2015
month.temp <- 9
# modify your working directory with the data file if necessary (using the setwd() command):
load("data_full.RData")  # full data without missing
data_ori <- data_DF
train <- subset(data_ori, year == year.temp & month == month.temp)
loc <- cbind(train$lon, train$lat)

c.temp <- c <- 0.6
# gamma estimation 
gamma <- c()
for (i in 1:3503){ #i-th location
  #i <- 945
  obs.temp <- subset(data_ori, lon == loc[i, 1] & lat == loc[i, 2])
  sort.temp <- sort(obs.temp$BA, decreasing = TRUE)
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
gamma[309]
which.max(gamma)

loc[309,]

    
#pdf("gamma_BA.pdf", height = 20, width = 15)
map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
map("world", fill = T, col = "gray90", xlim = c(-125, -63), ylim = c(24, 51))
quilt.plot(loc, gamma, ny = 35, add = TRUE)
title("Extreme value index (BA)", cex.main = 2)

     
    
   