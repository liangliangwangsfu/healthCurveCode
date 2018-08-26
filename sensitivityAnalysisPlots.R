folder0 = "SensitivityAnalysis/LHIN15_186/realDataPatient_LHIN15_186_MIN-2_MAX2_fixedThredI-1_fixedThredII1_knotsOption1_distBetween2knots10_nBatch2000_nIterInOneBatch3"
load(file=paste(folder0,"/result.RData",sep=""))

round(newMax,2)

apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))
meanVec <- t(round(apply(result_beta[nstart:nend,],2,mean),2))
quantilesVec <- t(round(apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.975),na.rm=T)),2))
cbind(t(meanVec), quantilesVec)



nSel <- 4
pdf(paste(folder0, "/selectedCurves_raw.pdf",sep=""), width=6,height=4)
#oneSample <- sample(1:nrow(mean.mat),nSel)
oneSample <- 20+(1:nSel)
colvec <- 1:nSel#c(1,4,2,6)
ltyvec <- 1:nSel#c(2,3,1,4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
nLevels = 7
title <- NULL
plot(
  c(0, 100), c(0, nLevels-1), xlim = c(1,90),type = "n", yaxt = 'n', main =
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
dev.off()


pdf(paste(folder0, "/selectedCurves_Setting1.pdf",sep=""), width=6,height=4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
oneCurve <- (colMeans(hi_resultList_0[[oneSample[1]]])-MIN)/(newMax-MIN)
plot(
  oneCurve, type = 'l', main = "",xlab = "",xlim = c(1,m), ylim =
    c(0,1), col = colvec[1], lwd = 2,lty = ltyvec[1],ylab = ""
)
for (j in 2:length(oneSample))
{
  oneCurve <-  (colMeans(hi_resultList_0[[oneSample[j]]])-MIN)/(newMax-MIN)
  if(lastDays[oneSample[j]]!=m) 
    oneCurve <- c(oneCurve, rep(0,m-lastDays[oneSample[j]]))
  lines(
    oneCurve, type = 'l', col = colvec[j], lwd = 2,lty = ltyvec[j]
  )
}
dev.off()



folder1 = "SensitivityAnalysis/LHIN15_186/realDataPatient_LHIN15_186_MIN-5_MAX5_fixedThredI-5_fixedThredII1_knotsOption1_distBetween2knots10_nBatch2000_nIterInOneBatch3"
load(file=paste(folder1,"/result.RData",sep=""))
round(newMax,2)

apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))
meanVec <- t(round(apply(result_beta[nstart:nend,],2,mean),2))
quantilesVec <- t(round(apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.975),na.rm=T)),2))
cbind(t(meanVec), quantilesVec)


nSel <- 4
pdf(paste(folder1, "/selectedCurves_raw.pdf",sep=""), width=6,height=4)
#oneSample <- sample(1:nrow(mean.mat),nSel)
oneSample <- 20+(1:nSel)
colvec <- 1:nSel#c(1,4,2,6)
ltyvec <- 1:nSel#c(2,3,1,4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
nLevels = 7
title <- NULL
plot(
  c(0, 100), c(0, nLevels-1), xlim = c(1,90),type = "n", yaxt = 'n', main =
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
dev.off()


pdf(paste(folder1, "/selectedCurves_Setting2.pdf",sep=""), width=6,height=4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
oneCurve <- (colMeans(hi_resultList_0[[oneSample[1]]])-MIN)/(newMax-MIN)
plot(
  oneCurve, type = 'l', main = "",xlab = "",xlim = c(1,m), ylim =
    c(0,1), col = colvec[1], lwd = 2,lty = ltyvec[1],ylab = ""
)
for (j in 2:length(oneSample))
{
  oneCurve <-  (colMeans(hi_resultList_0[[oneSample[j]]])-MIN)/(newMax-MIN)
  if(lastDays[oneSample[j]]!=m) 
    oneCurve <- c(oneCurve, rep(0,m-lastDays[oneSample[j]]))
  lines(
    oneCurve, type = 'l', col = colvec[j], lwd = 2,lty = ltyvec[j]
  )
}
dev.off()



folder2 = "SensitivityAnalysis/LHIN15_186/realDataPatient_LHIN15_186_MIN-10_MAX10_fixedThredI-10_fixedThredII1_knotsOption1_distBetween2knots10_nBatch2000_nIterInOneBatch3"
load(file=paste(folder2,"/result.RData",sep=""))
round(newMax,2)

apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.5,0.975),na.rm=T))
meanVec <- t(round(apply(result_beta[nstart:nend,],2,mean),2))
quantilesVec <- t(round(apply(result_beta[nstart:nend,],2,function(x)quantile(x,c(0.025,0.975),na.rm=T)),2))
cbind(t(meanVec), quantilesVec)


nSel <- 4
pdf(paste(folder2, "/selectedCurves_raw.pdf",sep=""), width=6,height=4)
#oneSample <- sample(1:nrow(mean.mat),nSel)
oneSample <- 20+(1:nSel)
colvec <- 1:nSel#c(1,4,2,6)
ltyvec <- 1:nSel#c(2,3,1,4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
nLevels = 7
title <- NULL
plot(
  c(0, 100), c(0, nLevels-1), xlim = c(1,90),type = "n", yaxt = 'n', main =
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
dev.off()


pdf(paste(folder2, "/selectedCurves_Setting3.pdf",sep=""), width=6,height=4)
par(mar = c(2, 2, 0.5, 0.5), mfrow = c(1,1))
oneCurve <- (colMeans(hi_resultList_0[[oneSample[1]]])-MIN)/(newMax-MIN)
plot(
  oneCurve, type = 'l', main = "",xlab = "",xlim = c(1,m), ylim =
    c(0,1), col = colvec[1], lwd = 2,lty = ltyvec[1],ylab = ""
)
for (j in 2:length(oneSample))
{
  oneCurve <-  (colMeans(hi_resultList_0[[oneSample[j]]])-MIN)/(newMax-MIN)
  if(lastDays[oneSample[j]]!=m) 
    oneCurve <- c(oneCurve, rep(0,m-lastDays[oneSample[j]]))
  lines(
    oneCurve, type = 'l', col = colvec[j], lwd = 2,lty = ltyvec[j]
  )
}
dev.off()


