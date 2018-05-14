rm(list = ls())
seed= 117
set.seed(seed)  # 1 argument: the seed

folder=paste("./Result-", seed, "/",sep="")
#folder=paste("./",sep="")
library("fda")
source("functions.R")
source("MCMC.R")
library(numDeriv)
library(fda)
library(MASS)
library(truncnorm)
library(mvtnorm)
library(MCMCpack)

beta_true <- c(8, 0.2,-0.1)
m=50
n_simu = 100
xVec=seq(1,m,length.out = m)

meanCurveFun <- function(x,m)
{
  ((x-m*0.3)/(2*m))^2+0.4
}
meanCurve = meanCurveFun(xVec,m)
plot(xVec,meanCurve, type="l",ylim=c(0,1))

FPC1 <- function(x,m){
  0.1-((x-m*0.5))^2/(m*400) + x*0.001
}
fpc1 <- FPC1(xVec,m)
plot(xVec,fpc1, type="l")

FPC2 <- function(x){
  0.15*cos(x*3.14/m)
}
fpc2 <- FPC2(xVec)
plot(xVec,fpc2, type="l")


FPC3 <- function(x){
  0.15*cos(x*6.28/m)+0.05
}
fpc3 <- FPC3(xVec)
plot(xVec,fpc3, type="l")

scoreVar <- c(3, 0.5, 0.3)
scoreSd <- sqrt(scoreVar)

scoreMat <- matrix(0,n_simu,3)
for(i in 1:3) scoreMat[,i] <- rnorm(n_simu,0,scoreSd[i])

data_simu <- matrix(0,n_simu,m)
for(i in 1:n_simu)
  data_simu[i,] <- scoreMat[i,1]*fpc1+scoreMat[i,2]*fpc2+scoreMat[i,3]*fpc3
matplot(t(data_simu),type = "l", xlab="Time")

MAX=10
MIN=-10
#fixedThred2 = 1

TTransform <- function(x)
{
  (x-MIN)/(MAX-MIN)
}

inverseTTransform <- function(x)
{
  x*(MAX-MIN)+MIN
}

data_simu <-data_simu + matrix(rep(meanCurve,n_simu),n_simu,m,byrow=T)
matplot(t(data_simu),type = "l", xlab="Time")
data_simu[data_simu>1]=1
data_simu[data_simu<0]=0
matplot(t(data_simu),type = "l", xlab="Time")

true_gamma <- c(-1, -2, -0.4, -0.5)
kappa(true_gamma)

data_simu=inverseTTransform(data_simu)
matplot(t(data_simu),type = "l", xlab="Time")

mean_data_simu <- apply(data_simu,1,function(x) mean(x, na.rm=T))
data_simu = data_simu - matrix(rep(mean_data_simu,m), n_simu,m,byrow = FALSE)
matplot(t(data_simu),type = "l", xlab="Time")

true_sigma2_epsilon = 0.8
Vsqrt = sqrt(true_sigma2_epsilon)*diag(rep(1,m))
R=diag(rep(1,m))
sigma2_epsilon_Mat <- mvrnorm(n_simu, rep(0,m), Vsqrt%*%R%*%Vsqrt)

true_sigma2_b <- 4
true_bvec <- rnorm(n_simu,0, sqrt(true_sigma2_b))

x0 <- rep(1,n_simu)
x1 <- rep(1,n_simu);

female <- sample(1:n_simu, round(0.5 * n_simu), replace = FALSE)
x1[female] = 0
###simulate age##
x2 <- sample(35:100,n_simu,replace = T)
X_sim <- cbind(1,x1,x2)
###simulate spline functions###
ncovariates = ncol(X_sim)
####simulate z####
threshvals <- kappa(true_gamma)
Z_sim <- matrix(NA,nr = n_simu,nc = m)
health.mat <- matrix(NA,nr = n_simu,nc = m)
#  norder = 4
for (i in 1:n_simu) {
  health.mat[i,]<-apply(rbind(data_simu[i,]+true_bvec[i], rep(MIN,m)),2,max)
  zerosInd=which(health.mat[i,]<=MIN)
  if(length(zerosInd)>0 && zerosInd[1]<m)
    health.mat[i,(zerosInd[1]+1):m]=MIN
  Z_sim[i,] <-health.mat[i,]  + rep(X_sim[i,1:3]%*% beta_true[1:3],m) + sigma2_epsilon_Mat[i,]
}
par(mfrow=c(1,2))
matplot(t(health.mat),type = "l", xlab="Time")  
matplot(t(Z_sim),type = "l", xlab="Time")  

y_sim <- apply(Z_sim, 1:2, function(x)findCat(x,threshvals))
sum(y_sim==0)
y_sim[health.mat<=MIN]=0

firstDeath_obs=apply(y_sim,1,function(x) which(x==0)[1])
firstDeath_obs[is.na(firstDeath_obs)]=-1

firstDeath_health=apply(health.mat,1,function(x) which(x==MIN)[1])
firstDeath_health[is.na(firstDeath_health)]=-1

goodSimu=firstDeath_obs==firstDeath_health
y_sim=y_sim[goodSimu,]
health.mat=health.mat[goodSimu,]
X_sim=X_sim[goodSimu,]

for (i in 1:n_simu) {
  zerosInd=which(y_sim[i,]==0)
  if(length(zerosInd)>0 && zerosInd[1]<m)
    y_sim[i,(zerosInd[1]+1):m]=0
}
sum(y_sim==0)
health.mat[y_sim==0]=MIN

X=X_sim
y=y_sim

n=nrow(X)
ncovariates = ncol(X)
d = rep(-1,n)
for (j in 1:n) {
  deadDays = which(y[j,] == 0)
  if (length(deadDays) > 0)
    d[j] = min(deadDays)
}
lastDays <- d
lastDays[d == -1] = m

X=X[lastDays>2,]
y=y[lastDays>2,]
Z_sim= Z_sim[lastDays>2,]
health.mat=health.mat[lastDays>2,]

n=nrow(X)
m=ncol(y)

ncovariates = ncol(X)
d = rep(-1,n)
for (j in 1:n) {
  deadDays = which(y[j,] == 0)
  if (length(deadDays) > 0)
    d[j] = min(deadDays)
}
lastDays <- d
lastDays[d == -1] = m

source("spline_basis.R")

nIterInOneBatch=2
nBatch=1000

a_lambda = 1/100
b_lambda = 100

temp = rgamma(100000, shape =  a_lambda, scale = b_lambda)
mean(temp)  
var(temp)

# mean and variance of the prior (normal) of \beta
mu_0 = rep(0,ncovariates)
sigma_0 <- diag(rep(10,ncovariates))
sd_gamma_hat = c(0.01, 0.01, 0.01, 0.01)

MAX=20
MIN=-20
runmcmcre<-oneRunMCMC(X, y, folder, nIterInOneBatch, nBatch,a_lambda, b_lambda, sd_gamma_hat)

source("displayDataFuns.R")

filename_lambda = paste(folder,"lambda",sep="")
filename_lambdamu = paste(folder,"lambdamu",sep="")  
filename_beta = paste(folder,"beta",sep="") 
filename_c = paste(folder,"c",sep="")
filename_cmu = paste(folder,"cmu",sep="")  
filename_gamma = paste(folder,"gamma",sep="")  
filename_accept = paste(folder,"accept",sep="")   
filename_sigma2_epsilon = paste(folder,"sigma2_epsilon",sep="")  
filename_sigma2_b = paste(folder,"sigma2_b",sep="")
filename_bvec = paste(folder,"bvec",sep="")  


result_bvec <- read.table(file = paste(filename_bvec,".txt",sep=""), header=FALSE) 
nend=nrow(result_bvec)
nstart=round(nend*0.3)


pdf(paste(filename_bvec,".pdf",sep=""),width=4*3,height=4)
par(mfrow=c(2,2*3),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
for(i in 1:12)
{
  temp = c(result_bvec[nstart:nend,i],true_bvec[i])
  hist(result_bvec[nstart:nend,i],main='b', xlim=range(temp))
  abline(v=true_bvec[i], lwd=2, col=2)
}
dev.off()

selSample=sample(1:n, 8)
pdf(paste(filename_bvec,"_rand8.pdf",sep=""),width=12,height=8)
par(mfrow=c(2,4),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)

for(i in selSample)
{
  temp = c(result_bvec[nstart:nend,i],true_bvec[i])
  hist(result_bvec[nstart:nend,i],main="", xlim=range(temp),col=4, xlab="", ylab="")
  abline(v=true_bvec[i], lwd=2, col=2)
}
dev.off()

result_sigma2_b <- read.table(file =  paste(filename_sigma2_b,".txt",sep=""), header=FALSE);
pdf(paste(filename_sigma2_b,"_trace.pdf",sep=""),width=9,height=9)
plot(result_sigma2_b[nstart:nend, 1],ylab='',xlab='',type='l')
dev.off()
pdf(paste(filename_sigma2_b,"_hist.pdf",sep=""),width=9,height=9)
hist(result_sigma2_b[nstart:nend, 1],col=4, xlab="", ylab="",main="")
abline(v=true_sigma2_b, lwd=2, col=2)
dev.off()

result_lambda <- read.table(file =  paste(filename_lambda,".txt",sep=""), header=FALSE);
nrows_lambda <- nrow(result_lambda)
ncols_lambda <- ncol(result_lambda)
pdf(paste(filename_lambda,".pdf",sep=""),width=9,height=9)
plot(result_lambda[nstart:nend, 1],ylab='',xlab='',type='l')
dev.off()

result_sigma2_epsilon <- read.table(file =  paste(filename_sigma2_epsilon,".txt",sep=""), header=FALSE);
pdf(paste(filename_sigma2_epsilon,"_trace.pdf",sep=""),width=9,height=9)
plot(result_sigma2_epsilon[,1], ylab="sigma2_epsilon")
dev.off()

pdf(paste(filename_sigma2_epsilon,"_traceAfterConverge.pdf",sep=""),width=9,height=9)
plot(result_sigma2_epsilon[nstart:nend,1],ylab='',type='l',xlab='')
dev.off()

pdf(paste(filename_sigma2_epsilon,"_hist.pdf",sep=""),width=9,height=9)
hist(result_sigma2_epsilon[nstart:nend,1], main="",xlab="",  ylab="",col=4)
abline(v=0.8,lwd=2,col=2)
dev.off()

result_gamma <- read.table(file = paste(filename_gamma,".txt",sep=""), header=FALSE)
pdf(paste(filename_gamma,"_alltrace.pdf",sep=""),width=10,height=5)
par(mfrow=c(2,2))
for(i in 1:4)
  plot(result_gamma[,i],ylab='',type='l',xlab='')
dev.off()

pdf(paste(filename_gamma,"_traceAfterConverge.pdf",sep=""),width=10,height=5)
par(mfrow=c(2,2))
for(i in 1:4)
  plot(result_gamma[nstart:nend,i],ylab='',type='l',xlab='')
dev.off()

pdf(paste(filename_beta,"_trace.pdf",sep=""),width=16,height=16)
result_beta <- read.table(file = paste(filename_beta,".txt",sep=""), header=FALSE) 
par(mfrow=c(2,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
for(i in 1:3) 
{
  plot(result_beta[,i],ylab='beta')
  abline(h=beta_true[i], lwd=2, col=2)
}
dev.off()

pdf(paste(filename_beta,"_hist.pdf",sep=""),width=8,height=4)
result_beta <- read.table(file = paste(filename_beta,".txt",sep=""), header=FALSE) 
par(mfrow=c(1,3))
for(i in 1:3) 
{
  temp = c(result_beta[nstart:nend,i], beta_true[i])
  temp = temp[!is.na(temp)]
  hist(result_beta[nstart:nend,i],xlim=range(temp),col=4,xlab="",main="")
  abline(v=beta_true[i], lwd=2, col=2)
  mtext(side=1,line=3,bquote(beta[.(i-1)]),cex=1) 
}
dev.off()

colMeans(result_beta[nstart:nend,])
beta_true

###hist_gamma###
burnin=nstart-1
gname = paste(filename_gamma,"_hist.pdf",sep="")  
pdf(gname,width=8,height=4)
result_gamma <- read.table(file = paste(filename_gamma,".txt",sep=""), header=FALSE)
par(mfrow=c(2,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
hist(result_gamma[-(1:burnin),1],ylab='',xlab='',main ="", col=4)
abline(v=true_gamma[1],lwd=3,col=2)
mtext(side=1,line=2,expression(gamma[1]),cex=1.2) 
hist(result_gamma[-(1:burnin),2],ylab='',xlab='',main ="", col=4)
abline(v=true_gamma[2],lwd=3,col=2)
mtext(side=1,line=2,expression(gamma[2]),cex=1.2) 
hist(result_gamma[-(1:burnin),3],ylab='',xlab='',main ="", col=4)
abline(v=true_gamma[3],lwd=3,col=2)
mtext(side=1,line=2,expression(gamma[3]),cex=1.2) 
hist(result_gamma[-(1:burnin),4],ylab='',xlab='',main ="", col=4)
abline(v=true_gamma[4],lwd=3,col=2)
mtext(side=1,line=2,expression(gamma[4]),cex=1.2) 
dev.off()

c_i_list=vector("list",n)           
hi_resultList=vector("list",n)           
hi_resultList_low=vector("list",n)
hi_resultList_upper=vector("list",n)   
for(i in 1:n)
{
  c_i <- read.table(file = paste(filename_c,"_",i, ".txt", sep=""), header=FALSE);
  c_i_list[[i]] <- c_i[burnin:nend,1:nbasisVec[i]]-colMeans(c_i[burnin:nend,1:nbasisVec[i]])+result_bvec[burnin:nend,i]
  hi_resultList[[i]]<-basismatList[[i]]%*% colMeans(c_i_list[[i]])
  hi_resultList_low[[i]]<-basismatList[[i]]%*% apply(c_i_list[[i]],2,function(x) quantile(x,0.025))
  hi_resultList_upper[[i]]<-basismatList[[i]]%*% apply(c_i_list[[i]],2,function(x) quantile(x,0.975))
}

for(i in 1:n)
{
  pdf(paste(filename_c,"_",i, ".pdf", sep=""),width=10,height=2)
  par(mfrow=c(1,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,3,0.2,0.5),cex.axis=1,las=1,
      mgp=c(1.5,0.5,0),adj=0.5)
  plot((hi_resultList[[i]]-MIN)/(MAX-MIN), type='l', ylab=paste("Health (", i,")",sep=""),xlab="Days",xlim=c(1,m), ylim=c(0,1))
  hilower=(hi_resultList_low[[i]]-MIN)/(MAX-MIN)
  hiupper=(hi_resultList_upper[[i]]-MIN)/(MAX-MIN)
  xseq = 1:length(hilower)
  polygon(c(rev(xseq), xseq), c(rev(hilower), hiupper), col = 'grey85', border = NA)
  lines(hilower, lwd=3, col=3, lty=3)
  lines(hiupper, lwd=3, col=3, lty=3)
  lines((health.mat[i,]-MIN)/(MAX-MIN), lwd=2, col=2, lty=2)
  lines((hi_resultList[[i]]-MIN)/(MAX-MIN),lwd=2,col=1)
  showOnePatientLocation(y[i,],1)
  dev.off()
}



