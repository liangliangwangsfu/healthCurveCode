rm(list = ls())

library("fda")
library("MCMCpack")
source("functions.R")

folderName="onesample-Nov6-MCMC"
folder=paste("./",folderName,"/",sep="")

load("allXy.RData")
set.seed(300)
nsub = 5000
indx = sample(1:nrow(X), nsub)
X=X[indx,]
y=y[indx,]
data<-removeDeathOnFirstDay(X,y)
X=data$X
y=data$y
lastDays=data$lastDays

source("spline_basis.R")
source("MCMC.R")
nIterInOneBatch=3
nBatch=5000

# hyper parameters for the smoonthing parameter \lambda
a_lambda = 1/100
b_lambda = 100
temp = rgamma(100000, shape =  a_lambda, scale = b_lambda)
mean(temp)  
var(temp)

# mean and variance of the prior (normal) of \beta
ncovariates=ncol(X)
mu_0 = rep(0,ncovariates)
sigma_0 <- diag(rep(10,ncovariates))
sd_gamma_hat = c(0.01, 0.01, 0.01, 0.01)
MAX=10
MIN=-10   #death
runmcmcre<-oneRunMCMC(X,y, folder, nIterInOneBatch, nBatch,a_lambda, b_lambda, sd_gamma_hat)
#save.image(file=paste(folderName,".RData",sep=""))

