rm(list=ls())

library(psych)
library(dplyr)
library(spdep)
library(caret)
library(patchwork)
library(splines2)
library(rstan)
rstan_options(auto_write =TRUE)
options(mc.cores = parallel::detectCores())

# modify your working directory with the data file if necessary (using the setwd() command):

# read data:
load("/Users/jongminlee/Desktop/EVA2021/Data/data_train.RData")
data_odd_DF <- readRDS("/Users/jongminlee/Desktop/EVA2021/Data/data_odd_DF.RDS")
data_odd_true <- data_train_DF %>% filter(year%%2 == 1) 

# show severity thresholds:
u_ba # for BA
u_cnt # for CNT
# show weights used for the weighted RPS:
weights_ba
weights_cnt

##### Benchmark ##### 
# remove rows with data to predict:
train_DF = data_odd_DF[!is.na(data_odd_DF$CNT),]
test_DF = data_odd_DF[is.na(data_odd_DF$CNT),]
# remove BA (which contains more NA values):
train_DF = subset(train_DF, select = -BA)
# in test data, remove BA and the NA column of CNT:
test_DF = subset(test_DF, select = -c(CNT, BA))
# train the model (Poisson regression with all covariates):
fit = glm(CNT ~ ., data = train_DF, family = poisson(link = "log"))
summary(fit)

# calculate estimated Poisson intensities (=means of predictive distributions):
pred_mean_cnt = predict(fit, test_DF, type = "response")
# calculate the matrix with estimated exceedance probability of the severity thresholds:
prediction_cnt = matrix(nrow = 41196, ncol = length(u_cnt))
for(k in 1:length(u_cnt)){
  prediction_cnt[,k] = ppois(u_cnt[k], lambda = pred_mean_cnt)
}
# prediction_cnt has to be submitted for the competition


# BA (log-normal regression): ####
# Procedure: we here fit the log-normal model only to positive values of BA, and then combine this with the count model above for the probability of BA=0.
# remove rows with data to predict:
train_DF = data_odd_DF[!is.na(data_odd_DF$BA),]
test_DF = data_odd_DF[is.na(data_odd_DF$BA),]
# remove CNT (which contains more NA values):
train_DF = subset(train_DF, select = -CNT)
# # in test data, remove CNT and the NA column of BA:
test_DF = subset(test_DF, select = -c(CNT, BA))
# use only positive BA values for fitting the log-normal GLM:
train_DF = train_DF[train_DF$BA > 0 ,]
# train the model:
fit = lm(log(BA) ~ ., data = train_DF) 
summary(fit)

# extract the predictive standard deviation:
sd_gauss = sd(fit$residuals)
sd_gauss

# extract predicted Gaussian means of BA>0:
pred_mean_logba = predict(fit, test_DF)
hist(pred_mean_logba)
# calculate estimated exceedance probabilities for the BA component:
prediction_ba = matrix(nrow = 41196, ncol = length(u_ba))
for(k in 1:length(u_ba)){
  # we here use prediction_cnt[,1] for the probability P(BA = 0).
  prediction_ba[,k] = prediction_cnt[,1] + (1-prediction_cnt[,1]) * pnorm(log(u_ba[k]), mean = pred_mean_logba, sd = sd_gauss)
}
# prediction_ba has to be submitted for the competition

##### score evaluation #####
# cnt 
cnt.score.ftn <- function(true.cnt, pred.cnt){
  true.cdf <- as.numeric(true.cnt<=u_cnt)
  est.cdf <- pred.cnt 
  return(sum(weights_cnt*(true.cdf-est.cdf)^2))
}
cnt.score <- vector(length=41196)
true.cnt = data_odd_true$CNT[is.na(data_odd_DF$CNT)]
for(i in 1:41196){
  cnt.score[i] <- cnt.score.ftn(true.cnt[i], prediction_cnt[i,])
}

# ba
ba.score.ftn <- function(true.ba, pred.ba){
  true.cdf <- as.numeric(true.ba<=u_ba)
  est.cdf <- pred.ba 
  return(sum(weights_ba*(true.cdf-est.cdf)^2))
}
ba.score <- vector(length=41196)
true.ba = data_odd_true$BA[is.na(data_odd_DF$BA)]

for(i in 1:41196){
  ba.score[i] <- ba.score.ftn(true.ba[i], prediction_ba[i,])
}

total.score <- sum(cnt.score + ba.score)
#benchmark score is : cnt = 2286.977, ba = 1475.452
