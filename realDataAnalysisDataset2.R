
#uniqueLevel = unique(X[,"Patient_LHIN"])
uniqueLevel <- c(2, 10,  9,  4,  1,  8,  6,  7,  11,  3,  12,  13,  15,   5,  14,  16)
nObsEachLevel <- matrix(c(2889, 1592, 4410, 4757, 2073, 4349, 2900, 3631, 3110, 2222, 1200, 1762, 186, 1921, 970,  4),1)
#nObsEachLevel <- matrix(rep(0,length(uniqueLevel)), 1)
colnames(nObsEachLevel) <- uniqueLevel 
#for(i in 1:length(uniqueLevel))
#  nObsEachLevel[i] <- sum(X[,"Patient_LHIN"]==uniqueLevel[i])

#for(k in 1: length(uniqueLevel)-1)
k=13
#while(k<length(uniqueLevel))
{
removelist <- ls()
removelist <- removelist[removelist!="uniqueLevel" || removelist!="nObsEachLevel" || removelist!="k"]
rm(removelist)
library("fda")
library("MCMCpack")
source("functions.R")

load("strokeData.RData")

folder=paste("./Results/realDataPatient_LHIN", uniqueLevel[k],"_", nObsEachLevel[k], "/",sep="")

#set.seed(200)
#nsub = 100
indx = X[,"Patient_LHIN"]==uniqueLevel[k]
nsub=sum(indx)
X=X[indx,1:2]
y=y[indx,]
data<-removeDeathOnFirstDay(X,y)
X=data$X
y=data$y
lastDays=data$lastDays

source("spline_basis.R")
source("MCMC.R")
nIterInOneBatch=5
nBatch=1000

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
#sd_gamma_hat = c(0.1, 0.1, 0.1, 0.1)
MAX=10
MIN=-10   #death
runmcmcre<-oneRunMCMC(X,y, folder, nIterInOneBatch, nBatch,a_lambda, b_lambda, sd_gamma_hat)
#save.image(file=paste(folderName,".RData",sep=""))
source("makeGraph.R")
k=k+1
}
