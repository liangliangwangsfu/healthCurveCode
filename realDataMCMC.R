folder=paste("./Results/realDataPatient_LHIN", uniqueLevel,"_", nObsThisLevel,"_MIN-", abs(MIN),"_MAX",MAX,"_nBatch",nBatch, "/",sep="")
indx = X[,"Patient_LHIN"]==uniqueLevel
nsub=sum(indx)
X=X[indx,1:2]
y=y[indx,]

X=X[1:50,]
y=y[1:50,]

data<-removeDeathOnFirstDay(X,y)
X=data$X
y=data$y
lastDays=data$lastDays
source("spline_basis.R")
#nIterInOneBatch=5
#nBatch=1000
# hyper parameters for the smoonthing parameter \lambda
a_lambda = 1/100
b_lambda = 100
temp = rgamma(100000, shape =  a_lambda, scale = b_lambda)
#mean(temp)  
#var(temp)
# mean and variance of the prior (normal) of \beta
ncovariates=ncol(X)
mu_0 = rep(0,ncovariates)
sigma_0 <- diag(rep(10,ncovariates))
sd_gamma_hat = c(0.01, 0.01, 0.01, 0.01)
#sd_gamma_hat = c(0.1, 0.1, 0.1, 0.1)
#MAX=10
#MIN=-10   #death
ptm=proc.time()
runmcmcre<-oneRunMCMC(X,y, folder, nIterInOneBatch, nBatch,a_lambda, b_lambda, sd_gamma_hat)
mcmcTime=proc.time()-ptm
source("makeGraph.R")
save.image(file=paste(folder,"/result.RData",sep=""))
