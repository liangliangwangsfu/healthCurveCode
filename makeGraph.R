source("displayDataFuns.R")
library("fda")

#folder="./onesample-May13-MCMC/"  
filename_lambda = paste(folder,"lambda",sep="")
filename_beta = paste(folder,"beta",sep="") 
filename_c = paste(folder,"c",sep="")
filename_cmu = paste(folder,"cmu",sep="")  
filename_gamma = paste(folder,"gamma",sep="")  
filename_accept = paste(folder,"accept",sep="")   
filename_sigma2_epsilon = paste(folder,"sigma2_epsilon",sep="")  
filename_sigma2_b = paste(folder,"sigma2_b",sep="")
filename_bvec = paste(folder,"bvec",sep="")  
result_bvec <- read.table(file = paste(filename_bvec,".txt",sep=""), header=FALSE) 
result_sigma2_b <- read.table(file =  paste(filename_sigma2_b,".txt",sep=""), header=FALSE);
result_lambda <- read.table(file =  paste(filename_lambda,".txt",sep=""), header=FALSE);
nrows_lambda <- nrow(result_lambda)
ncols_lambda <- ncol(result_lambda)
nend <- nrows_lambda
nstart <- round(nend*0.3)
burnin <- nstart
result_sigma2_epsilon <- read.table(file =  paste(filename_sigma2_epsilon,".txt",sep=""), header=FALSE);
result_gamma <- read.table(file = paste(filename_gamma,".txt",sep=""), header=FALSE)

c_i_list=vector("list",n)           
hi_resultList=vector("list",n)           
for(i in 1:n)
{
  c_i_list[[i]] <- read.table(file = paste(filename_c,"_",i, ".txt", sep=""), header=FALSE);
}

hi_resultList=vector("list",n)          
hi_resultList_lower=vector("list",n)
hi_resultList_upper=vector("list",n)          
for(i in 1:n)
{
  hi_resultList[[i]]<-basismatList[[i]]%*% colMeans(c_i_list[[i]][burnin:nend,1:nbasisVec[i]])
  hi_resultList_lower[[i]] <- basismatList[[i]]%*%apply(c_i_list[[i]],2,function(x)quantile(x,0.025)) 
  hi_resultList_upper[[i]] <- basismatList[[i]]%*%apply(c_i_list[[i]],2,function(x)quantile(x,0.975)) 
}

mean.mat = matrix(NA,n,m)
hi_resultList_0=vector("list",n)
for(i in 1:n)
{
  re.mat = matrix(NA,nend-burnin,m)
  for(j in 1:(nend-burnin))
  {
    re.mat[j,1:lastDays[i]] <- basismatList[[i]]%*% t(c_i_list[[i]][j+burnin,1:nbasisVec[i]])
    re.mat[j,1:lastDays[i]] <- re.mat[j,1:lastDays[i]] - mean(re.mat[j,1:lastDays[i]]) +result_bvec[j+burnin,i]
   #  re.mat[j,1:lastDays[i]] <- re.mat[j,1:lastDays[i]]  +result_bvec[j+burnin,i]
  }
  re.mat[re.mat<MIN]=MIN
  hi_resultList_0[[i]]=re.mat
  mean.mat[i,]=colMeans(re.mat,na.rm=T)
}

mean.mat.2 = matrix(NA,nend-burnin,m)
for(j in 1:(nend-burnin))
{
  re.mat = matrix(NA,n,m)
  for(i in 1:n)
  {
    #cbind(basismatList[[i]]%*% t(c_i_list[[i]][j+burnin,1:nbasisVec[i]])+result_bvec[j+burnin,i], basismatList[[i]]%*% t(c_i_list[[i]][j+burnin,1:nbasisVec[i]])- mean(re.mat[i,1:lastDays[i]])+result_bvec[j+burnin,i])
    re.mat[i,1:lastDays[i]] <- basismatList[[i]]%*% t(c_i_list[[i]][j+burnin,1:nbasisVec[i]])
    cat(mean(re.mat[i,1:lastDays[i]]), result_bvec[j+burnin,i], "\n")
    #re.mat[i,1:lastDays[i]] <- re.mat[i,1:lastDays[i]] - mean(re.mat[i,1:lastDays[i]]) +result_bvec[j+burnin,i]
    re.mat[i,1:lastDays[i]] <- re.mat[i,1:lastDays[i]] +result_bvec[j+burnin,i]
  }
  re.mat[re.mat<MIN]=MIN
  mean.mat.2[j,]=colMeans(re.mat,na.rm=T)
}


mean(unlist(result_bvec[burnin:nend,]))

#########################################################################################
newMax = quantile(mean.mat,0.995,na.rm=T)  
#newMax = quantile(rowMeans(mean.mat),0.99,na.rm=T)
#newMax = max((mean.mat), na.rm=T)
himean  <- (colMeans(mean.mat.2,na.rm = T)-MIN)/(newMax-MIN)
mean(himean)
hilower <- (apply(mean.mat.2,2,function(x)  quantile(x, 0.025, na.rm = T))-MIN)/(newMax-MIN)
hiupper <- (apply(mean.mat.2,2,function(x)  quantile(x, 0.975, na.rm = T))-MIN)/(newMax-MIN)
pdf(paste(folder, "meanCurve.pdf", sep=""),width=6,height=6)
par(mfrow=c(1,1))
plot(himean, type='l', main=paste("Health Curve"),xlab="Days", ylab="Health",xlim=c(1,m), ylim=c(0.4,0.9))
lines(hilower,lty=2)
lines(hiupper,lty=2)
xseq = 1:m
polygon(c(rev(xseq), xseq), c(rev(hilower), hiupper), col = 'grey85', border = NA)
lines(himean, type='l', lwd=2)
dev.off() 

scaledhMat = (mean.mat-MIN)/(newMax-MIN)
normHMat = scaledhMat 
normHMat[normHMat<0]=0
normHMat[is.na(normHMat)]=0 
normHMat[normHMat>1]=1

########################################################################################

meanH <- apply(normHMat, 2, function(x)  mean(x[x!=0]))
medianH<-apply(normHMat,2, function(x)  median(x[x!=0]))

normHMat[normHMat>1]=1
normHMat[normHMat<0]=0
par(mfrow = c(1, 1))
pdf(paste(folder,"normalizedHealthCurveWithMean.pdf",sep=""), width=6,height=6)
matplot(t(normHMat),type = "l", xlab="Time")
lines(meanH,type="l", col=1, lwd=3)
dev.off()

##################################################
normHMat[normHMat>1]=1
normHMat[normHMat<0]=0
pdf(paste(folder, "selectedCurves.pdf",sep=""), width=6,height=6)
nSel <- 4
#oneSample <- sample(1:nrow(mean.mat),nSel)
#oneSample <- 20+(1:nSel)
oneSample <- c(21, 22, 33, 98)
colvec <- 1:nSel#c(1,4,2,6)
ltyvec <- 1:nSel#c(2,3,1,4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(2,1))
nLevels = 7
title <- NULL
plot(
  c(0, 100), c(0, nLevels + 1), xlim = c(1,90),type = "n", yaxt = 'n', main =
    "",xlab = '', ylab = "",col = 1, lwd = 2
)
axis(
  2, 0:(nLevels - 1), 0:(nLevels - 1), font = 1,col = 1, lwd = 1
)
for (j in 1:length(oneSample)) {
  x_i <- 1:length(y[oneSample[j],])
  re <- stepLines(x_i,y[oneSample[j],])
  lines(re$xstep, re$ystep, col = colvec[j], lwd = 2,lty = ltyvec[j])
}

oneCurve <- (colMeans(hi_resultList_0[[oneSample[1]]])-MIN)/(newMax-MIN)
plot(
  oneCurve, type = 'l', main = "",xlab = "",xlim = c(1,m), ylim =
    c(0,1), col = colvec[1], lwd = 2,lty = ltyvec[1],ylab = ""
)
for (j in 2:length(oneSample))
{
  oneCurve <-  (colMeans(hi_resultList_0[[oneSample[j]]])-MIN)/(newMax-MIN)
  if(lastDays[oneSample[j]]!=m) 
    oneCurve <- c(oneCurve[1:lastDays[oneSample[j]]], rep(0,m-lastDays[oneSample[j]]))
  lines(
    oneCurve, type = 'l', col = colvec[j], lwd = 2,lty = ltyvec[j]
  )
}
dev.off()

###############Functional PCA#########################################################################

timepts = seq(1,m,1)
norder = 4 ## cubic B-spline
nbasis = norder + length(timepts) - 2; 
## create cubic B-spline basis functions, class 'basisfd'
spline.basis = create.bspline.basis(rangeval = c(1,m),nbasis,norder,timepts)
basismat   = eval.basis(timepts, spline.basis);
D2Lfd <- int2Lfd(m = 2)
D2fdPar   <- fdPar(spline.basis, D2Lfd, 1e4)
health.fd <- Data2fd(y=t(normHMat), timepts, spline.basis)
PCAobjects1 = pca.fd(health.fd, nharm = 3,harmfdPar = D2fdPar)
PCAobjects1$varprop
tobs = timepts

########################################################################
pdf(paste(folder, "scaledfpcaResult.pdf",sep=""), width=9,height=3)
par(mfrow = c(1, 3))
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,1],type = "l",xlab = "Time",
  ylab = paste(
    "FPC 1 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[1]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,2],type = "l",xlab = "Time",
  ylab = paste(
    "FPC 2 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[2]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,3],type = "l",xlab = "Time",
  ylab = paste(
    "FPC 3 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[3]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
dev.off()

#pdf(paste(folder,"PCscore.pdf",sep=""), width=6,height=6)
#par(mfrow=c(1,1))
#par(mar = c(4.5, 4.5, 0.5, 0.5), mfrow = c(1,1))
#pcscore1<-PCAobjects1$scores[,1]
#pcscore2<--PCAobjects1$scores[,2]
#plot(pcscore1[d==-1], pcscore2[d==-1], xlim=range(pcscore1), ylim=range(pcscore2), xlab="PC Score 1", ylab="PC Score 2")
#points(pcscore1[d>-1], pcscore2[d>-1], col=2, pch = 8)
#dev.off()

########################################################################################
pdf(paste(folder, "meanCurveScaledfpcaResult.pdf",sep=""), width=9,height=9)
par(mfrow = c(2, 2))
plot(himean, type='l', main=paste("Mean Health Curve"),xlab="Time (Days)", ylab="Health",xlim=c(1,m), ylim=c(0.7,0.9))
lines(hilower,lty=2)
lines(hiupper,lty=2)
xseq = 1:m
polygon(c(rev(xseq), xseq), c(rev(hilower), hiupper), col = 'grey85', border = NA)
lines(himean, type='l', lwd=2)
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,1],type = "l", lwd=2, xlab = "Time (Days)",
  ylab ="", main= paste(
    "FPC 1 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[1]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,2],type = "l", lwd=2,  xlab = "Time (Days)",
  ylab ="",  main=paste(
    "FPC 2 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[2]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
plot(
  tobs,basismat %*% PCAobjects1$harmonics$coefs[,3],type = "l", lwd=2,  xlab = "Time (Days)",
  ylab = "", main = paste(
    "FPC 3 (", sprintf("%.1f%%", 100 * PCAobjects1$varprop[3]),")",sep =
      ""
  )
)
lines(tobs,rep(0,length(tobs)),type = "l",lty = 2)
dev.off()

apply(PCAobjects1$scores,2,var)

fpc1 = basismat %*% PCAobjects1$harmonics$coefs[,1]
fpc2 = basismat %*% PCAobjects1$harmonics$coefs[,2]
fpc3 = basismat %*% PCAobjects1$harmonics$coefs[,3]

which(fpc1==max(fpc1))
which(fpc3>0)
# save(PCAobjects1,basismat,fpc1,fpc2,fpc3, file=paste(folder,"pcaResult.RData",sep=""))

##########################################################################################

#pdf(paste(folder,"PCscore2.pdf",sep=""), width=6,height=6)
#par(mfrow=c(1,1))
#par(mar = c(4.5, 4.5, 0.5, 0.5), mfrow = c(1,1))
#pcscore1<-PCAobjects1$scores[,1]
#pcscore2<--PCAobjects1$scores[,2]
#plot(pcscore1[d==-1], pcscore2[d==-1], xlim=range(pcscore1), ylim=range(pcscore2), xlab="PC Score 1", ylab="PC Score 2")
#points(pcscore1[d>-1], pcscore2[d>-1], col=2, pch = 8)
#dev.off()

############### K means clustering based on the first two PC scores ##################################################
#algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")
re.kmeans <- kmeans(PCAobjects1$scores[,1:2], 4, iter.max = 200, algorithm="MacQueen")
re <- re.kmeans$cluster
save(re.kmeans, file= paste(folder,"re.kmeans.RData",sep=""))

pdf(paste(folder,"PCscore1score2.pdf",sep=""), width=6,height=6)
pcscore1<-PCAobjects1$scores[,1]
pcscore2<--PCAobjects1$scores[,2]
plot(pcscore1[re==1], pcscore2[re==1], xlim=range(pcscore1), ylim=range(pcscore2), xlab="PC Score 1", ylab="PC Score 2")
points(pcscore1[re==2], pcscore2[re==2], col=2, pch = 2)
points(pcscore1[re==3], pcscore2[re==3], col=3, pch = 3)
points(pcscore1[re==4], pcscore2[re==4], col=4, pch = 4)
dev.off()

sum(re==1)/n  #[1] 0.05594973
sum(re==2)/n  #[1] 0.2150821
sum(re==3)/n  #[1] 0.5041557
sum(re==4)/n  #[1] 0.2248125

pdf(paste(folder,"clusters.pdf",sep=""), width=8,height=8)
par(mar = c(4.5, 4.5, 0.5, 0.5), mfrow = c(2,2))
quantile.thre = 0.05
dist1.vec = apply(PCAobjects1$scores[,1:2], 1, function(x) dist(rbind(x,re.kmeans$centers[1,])))
sel.ind = which(dist1.vec<quantile(dist1.vec, quantile.thre))
matplot(t(normHMat[sel.ind,]),type = "l", xlab="Time (Days)", col=1, ylim=c(0,1), lty=1)
lines(apply(normHMat[sel.ind,], 2, mean), lty = 4, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.975))), lty = 2, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.025))), lty = 2, col = 6, lwd = 3)


dist1.vec = apply(PCAobjects1$scores[,1:2], 1, function(x) dist(rbind(x,re.kmeans$centers[2,])))
sel.ind = which(dist1.vec<quantile(dist1.vec, quantile.thre))
matplot(t(normHMat[sel.ind,]),type = "l", xlab="Time (Days)", col=2, ylim=c(0,1), lty=1)
lines(apply(normHMat[sel.ind,], 2, mean), lty = 4, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.975))), lty = 2, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.025))), lty = 2, col = 6, lwd = 3)



dist1.vec = apply(PCAobjects1$scores[,1:2], 1, function(x) dist(rbind(x,re.kmeans$centers[3,])))
sel.ind = which(dist1.vec<quantile(dist1.vec, quantile.thre))
matplot(t(normHMat[sel.ind,]),type = "l", xlab="Time (Days)", col=3, ylim=c(0,1), lty=1)
lines(apply(normHMat[sel.ind,], 2, mean), lty = 4, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.975))), lty = 2, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.025))), lty = 2, col = 6, lwd = 3)


dist1.vec = apply(PCAobjects1$scores[,1:2], 1, function(x) dist(rbind(x,re.kmeans$centers[4,])))
sel.ind = which(dist1.vec<quantile(dist1.vec, quantile.thre))
matplot(t(normHMat[sel.ind,]),type = "l", xlab="Time (Days)", col=4, ylim=c(0,1), lty=1)
lines(apply(normHMat[sel.ind,], 2, mean), lty = 4, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.975))), lty = 2, col = 6, lwd = 3)
lines(apply(normHMat[sel.ind,], 2, function(x)(quantile(x, probs = 0.025))), lty = 2, col = 6, lwd = 3)


dev.off()

##########################################################################################

#pdf(paste(folder,"PCscore.pdf",sep=""), width=6,height=6)
#par(mfrow=c(1,1))
#par(mar = c(4.5, 4.5, 0.5, 0.5), mfrow = c(1,1))
#pcscore1<-PCAobjects1$scores[,1]
#pcscore2<--PCAobjects1$scores[,2]
#plot(pcscore1[d==-1], pcscore2[d==-1], xlim=range(pcscore1), ylim=range(pcscore2), xlab="PC Score 1", ylab="PC Score 2")
#points(pcscore1[d>-1], pcscore2[d>-1], col=2, pch = 8)
#dev.off()

##########################################################################################

pdf(paste(filename_sigma2_b,"_trace.pdf",sep=""),width=9,height=9)
plot(result_sigma2_b[nstart:nend, 1],ylab='',xlab='',type='l')
dev.off()
pdf(paste(filename_sigma2_b,"_hist.pdf",sep=""),width=9,height=9)
hist(result_sigma2_b[nstart:nend, 1],ylab='',xlab='',main ="", col=4)
dev.off()
round(c(mean(result_sigma2_b[nstart:nend, 1]), quantile(result_sigma2_b[nstart:nend, 1], c(0.025, 0.975))),3)
#24.034 23.108 25.019 

pdf(paste(filename_lambda,".pdf",sep=""),width=9,height=9)
plot(result_lambda[nstart:nend, 1],ylab='',xlab='',type='l')
dev.off()

pdf(paste(filename_lambda,"_hist.pdf",sep=""),width=9,height=9)
hist(result_lambda[nstart:nend, 1],ylab='',xlab='',main ="", col=4)
dev.off()
round(c(mean(result_lambda[nstart:nend, 1]), quantile(result_lambda[nstart:nend, 1], c(0.025, 0.975))),3)
#1.114 1.089 1.139 

pdf(paste(filename_sigma2_epsilon,"_trace.pdf",sep=""),width=9,height=9)
plot(result_sigma2_epsilon[,1], ylab="sigma2_epsilon")
dev.off()

pdf(paste(filename_sigma2_epsilon,"_traceAfterConverge.pdf",sep=""),width=9,height=9)
plot(result_sigma2_epsilon[nstart:nend,1],ylab='',type='l',xlab='')
dev.off()

pdf(paste(filename_sigma2_epsilon,"_hist.pdf",sep=""),width=9,height=9)
hist(result_sigma2_epsilon[nstart:nend,1],ylab='',xlab='',main ="", col=4)
dev.off()

round(c(mean(result_sigma2_epsilon[nstart:nend,1]), quantile(result_sigma2_epsilon[nstart:nend,1], c(0.025, 0.975))),3)
#0.7426141 0.7267414 0.7590022 

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
par(mfrow=c(3,1),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
for(i in 1:3) 
{
  plot(result_beta[,i],ylab='beta')
}
dev.off()

pdf(paste(filename_beta,"_hist.pdf",sep=""),width=16,height=8)
result_beta <- read.table(file = paste(filename_beta,".txt",sep=""), header=FALSE) 
par(mfrow=c(1,3),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
for(i in 1:3) 
{
  temp = result_beta[nstart:nend,i]
  temp = temp[!is.na(temp)]
  hist(result_beta[nstart:nend,i],xlim=range(temp),ylab='',xlab='',main ="", col=4)
}
dev.off()


apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))
meanVec <- t(round(apply(result_beta[nstart:nend,],2,mean),3))
quantilesVec <- t(round(apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.975),na.rm=T)),3))
cbind(t(meanVec), quantilesVec)
#      2.5%  97.5%
#  V1  3.234  3.101  3.367
#  V2  0.114  0.074  0.156
#  V3 -0.061 -0.063 -0.059


pdf(paste(filename_beta,"_hist_1.pdf",sep=""),width=4,height=4)
result_beta <- read.table(file = paste(filename_beta,".txt",sep=""), header=FALSE) 
par(mfrow=c(1,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
i=2
  temp = result_beta[nstart:nend,i]
  temp = temp[!is.na(temp)]
  hist(result_beta[nstart:nend,i],xlim=range(temp),ylab='',xlab='',main ="", col=4)
dev.off()


pdf(paste(filename_beta,"_hist_2.pdf",sep=""),width=4,height=4)
result_beta <- read.table(file = paste(filename_beta,".txt",sep=""), header=FALSE) 
par(mfrow=c(1,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
i=3
temp = result_beta[nstart:nend,i]
temp = temp[!is.na(temp)]
hist(result_beta[nstart:nend,i],xlim=range(temp),ylab='',xlab='',main ="", col=4)
dev.off()


###hist_gamma###
gname = paste(filename_gamma,"_hist.pdf",sep="")  
pdf(gname,width=8,height=4)
result_gamma <- read.table(file = paste(filename_gamma,".txt",sep=""), header=FALSE)
par(mfrow=c(2,2),oma=c(0.2,1.75,0.2,0.2),mar=c(3.5,2,0.2,3),cex.axis=1,las=1,
    mgp=c(1,0.5,0),adj=0.5)
hist(result_gamma[nstart:nend,1],ylab='',xlab='',main ="", col=4)
mtext(side=1,line=2,expression(gamma[1]),cex=1.2) 
hist(result_gamma[nstart:nend,2],ylab='',xlab='',main ="", col=4)
mtext(side=1,line=2,expression(gamma[2]),cex=1.2) 
hist(result_gamma[nstart:nend,3],ylab='',xlab='',main ="", col=4)
mtext(side=1,line=2,expression(gamma[3]),cex=1.2) 
hist(result_gamma[nstart:nend,4],ylab='',xlab='',main ="", col=4)
mtext(side=1,line=2,expression(gamma[4]),cex=1.2) 
dev.off()

t(round(apply(result_gamma[nstart:nend,],2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T)),3))
# 2.5%    50%  97.5%
# V1 -0.046 -0.025 -0.006
# V2 -2.730 -2.692 -2.657
# V3 -0.757 -0.745 -0.733
# V4 -0.576 -0.562 -0.547

#save.image(file= paste(folder,"result.RData",sep=""))

