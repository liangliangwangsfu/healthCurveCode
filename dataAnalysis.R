#! /usr/bin/env Rscript

# Input arguments:
#./dataAnalysis.R  "k=1" "nIterInOneBatch=2"   "nBatch=20" "nIterInOneBatch=3" "MAX=100" "MIN=-5" "fixedThred1 = -1" "fixedThred2 = 1" "knotsOption = 1" "distBetween2knots = 10"

rm(list=ls())
library("fda")
library("MCMCpack")
library(numDeriv)
library(fda)
library(MASS)
library(truncnorm)
source("functions.R")
source("MCMC.R")

args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text=args[i]))
}

uniqueLevel = c(2, 10,  9,  4,  1,  8,  6,  7,  11,  3,  12,  13,  15,   5,  14,  16)
nObsEachLevel = matrix(c(2889, 1592, 4410, 4757, 2073, 4349, 2900, 3631, 3110, 2222, 1200, 1762, 186, 1921, 970,  4),1)
colnames(nObsEachLevel) = uniqueLevel 
load("strokeData.RData")
uniqueLevel=uniqueLevel[k]
nObsThisLevel=nObsEachLevel[k]

#Create a folder to store the results and plots
folder = paste("./Results/realDataPatient_LHIN", uniqueLevel,"_", nObsThisLevel,"_MIN-", abs(MIN),"_MAX",MAX,"_fixedThredI-",abs(fixedThred1), "_fixedThredII",fixedThred2,  "_knotsOption", knotsOption,"_distBetween2knots",distBetween2knots, "_nBatch",nBatch, "_nIterInOneBatch",nIterInOneBatch, "/",sep="")
if (!file.exists(folder)){
  dir.create(folder)
}
logFile = paste(folder, "/log.txt",sep="")
cat("\n Patient_LHIN", uniqueLevel, "\n nObs", nObsThisLevel,"\n MIN", MIN,"\n MAX",MAX,"\n fixedThredI",fixedThred1, "\n fixedThredII",fixedThred2,  "\n knotsOption", knotsOption,"\n distBetween2knots",distBetween2knots, "\n nBatch",nBatch, "\n nIterInOneBatch",nIterInOneBatch, file = logFile, sep = " ", append = FALSE)

file.copy(from="dataAnalysis.R",to=paste(folder, "/dataAnalysis.R",sep=""))
indx = X[,"Patient_LHIN"]==uniqueLevel
nsub=sum(indx)
X=X[indx,1:2]
y=y[indx,]

# set.seed(300)
# nsub = 50
# indx = sample(1:nrow(X), nsub)
# X=X[indx,1:2]
# y=y[indx,]

data<-removeDeathOnFirstDay(X,y)
X=data$X
y=data$y
lastDays=data$lastDays
source("spline_basis.R")
# hyper parameters for the smoonthing parameter \lambda
a_lambda = 1/100
b_lambda = 100
temp = rgamma(100000, shape =  a_lambda, scale = b_lambda)
#mean(temp)  
#var(temp)
# mean and variance of the prior (normal) of \beta
ncovariates=ncol(X)
mu_0 = rep(0,ncovariates)
sigma_0 = diag(rep(10,ncovariates))
sd_gamma_hat = c(0.01, 0.01, 0.01, 0.01)

ptm=proc.time()
cat("\n start time: ", ptm, file = logFile, sep = " ", append = TRUE)
runmcmcre=oneRunMCMC(X,y, folder, nIterInOneBatch, nBatch,a_lambda, b_lambda, sd_gamma_hat)
endtime <- proc.time()
cat("\n end time: ", endtime, file = logFile, sep = " ", append = TRUE)
mcmcTime=endtime-ptm
cat("\n total MCMC time: ", mcmcTime, file = logFile, sep = " ", append = TRUE)
source("makeGraph.R")
save.image(file=paste(folder,"/result.RData",sep=""))
