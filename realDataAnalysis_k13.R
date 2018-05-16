rm(list=ls())
library("fda")
library("MCMCpack")
library(numDeriv)
library(fda)
library(MASS)
library(truncnorm)
source("functions.R")
source("MCMC.R")
#uniqueLevel = unique(X[,"Patient_LHIN"])
uniqueLevel <- c(2, 10,  9,  4,  1,  8,  6,  7,  11,  3,  12,  13,  15,   5,  14,  16)
nObsEachLevel <- matrix(c(2889, 1592, 4410, 4757, 2073, 4349, 2900, 3631, 3110, 2222, 1200, 1762, 186, 1921, 970,  4),1)
#nObsEachLevel <- matrix(rep(0,length(uniqueLevel)), 1)
colnames(nObsEachLevel) <- uniqueLevel 
#for(i in 1:length(uniqueLevel))
#  nObsEachLevel[i] <- sum(X[,"Patient_LHIN"]==uniqueLevel[i])
k=13
nIterInOneBatch=3
nBatch=100
MAX=10
MIN=-10
load("strokeData.RData")
#realDataMCMC(uniqueLevel[k], nObsEachLevel[k], nIterInOneBatch, nBatch, MAX, MIN)  
uniqueLevel=uniqueLevel[k]
nObsThisLevel=nObsEachLevel[k]
ptm=proc.time()
source("realDataMCMC.R")  
proc.time()-ptm
