
#setwd("/Users/shijiaw/Desktop/health_curve/RCode")

## Functions about data processing

stepLines <- function(x,y)
{
len.y <- length(y)
changePos <- (2:len.y)[y[2:len.y]!=y[1:(len.y-1)]]
if(length(changePos)==0)
return(list(xstep=x,ystep=y))
else
{
xnew <- insertMultiVals(x,x[changePos],changePos)
ynew <- insertMultiVals(y,y[changePos-1],changePos)
return(list(xstep=xnew,ystep=ynew))
}
}


jumpType <- function(x,y)
{
len.y <- length(y)
changePos <- (2:len.y)[y[2:len.y]!=y[1:(len.y-1)]]
changeNums <-length(changePos)
jumpMat <- NULL
if(changeNums>0)
{
jumpMat <- matrix(0,changeNums,2)
for(i in 1:length(changePos))
{
jumpMat[i,1] <- y[changePos[i]-1]
jumpMat[i,2] <- y[changePos[i]]
}
}
return(jumpMat)
}



insert <- function(v,e,pos){
  return(c(v[1:(pos-1)],e,v[(pos):length(v)]))
}


insertMultiVals <- function(v,iv,posv){
currentv <- v
for(i in length(posv):1)
{
currentv <- insert(currentv, iv[i], posv[i])
}
return(currentv) 
}


showPatientLocations <- function()
{
pt_vars <- read.csv("../../healthcurve/stroke_patients/pt_vars.csv")
stroke_location <- read.csv("../../healthcurve/stroke_patients/stroke_location.csv")
dim(stroke_location)

ptKeys <- unique(stroke_location$pat_cihi_key) 
npt <- length(ptKeys)  # patient counts

locationLevels <- levels(stroke_location[,3])
locationLevels[1]
nLevels <- length(locationLevels)
nLevels = 7
 

dev.new(width=10, height=4)
npt=50
for(i in 1:npt)
{  
 data_i <-stroke_location[stroke_location$pat_cihi_key==ptKeys[i],]
 ptInfo <- pt_vars[pt_vars$pat_cihi_key==ptKeys[i], 2:4]
 title <- paste(i, "/", npt,  "Region:", ptInfo[1], " Sex:", ptInfo[2], " Age:", ptInfo[3]) 
 par(las=1)
 par(mar=c(4.1,13.1, 4.1,4.1))
 plot(c(0, 100), c(0, nLevels+1), xlim=c(1,m),type = "n", yaxt='n', main=title, xlab='Days', ylab="")
 axis(2, 1:(nLevels+1), c(locationLevels, ""), font=1,col='black', lwd=1)
 x <- data_i[,2]
 y <- data_i[,3]
 re <- stepLines(x,y)
 lines(re$xstep, re$ystep)
}
}



showOnePatientLocation <- function(location, lineColor)
{ 
 nLevels=7
# par(mar=c(4.1,13.1, 4.1,4.1))
 #title <- paste("Patient  Age ", age) 
 title <- NULL
 plot(c(0, 100), c(0, nLevels+1), xlim=c(1,m),type = "n", yaxt='n', main=title,xlab='', ylab="",col=lineColor)
 axis(2, 0:(nLevels-1), 0:(nLevels-1), font=1,col=lineColor, lwd=1)
 x_i <- 1:length(location)
 re <- stepLines(x_i,location)
 lines(re$xstep, re$ystep)
 }






drawTransit <- function(){
    jumpMat <- matrix(0,nLevels,nLevels)
    for(i in 1:npt)
    {
        cat(i, " ")
        # data_i <- stroke_location[stroke_location$pat_cihi_key==ptKeys[i],]
        startIx <- (i-1)*90+1
        endIx <- i*90
        data_i <- stroke_location[startIx:endIx,]
        x <- data_i[,2]
        # print(x)
        y <- data_i[,3]
        jumps <- jumpType(x,y)
        jumpNum <- nrow(jumps)
        if(!is.null(jumpNum)){
            for(j in 1:jumpNum)
            jumpMat[jumps[j,1], jumps[j,2]] <-  jumpMat[jumps[j,1], jumps[j,2]]+1
        }
    }

jumpRows <- rowSums(jumpMat)
jumpPercent <- jumpMat/rowSums(jumpMat) #equal to for(i in 1:7) print(jumpMat[i,]/jumpRows[i])
transitMat <- round(jumpPercent, 2)
locationNames <- c("Acute","Dead","Emergency","Home","Chronic","Rehabilitation","Long-term")

rownames(transitMat) <- locationNames
colnames(transitMat) <- locationNames
library(diagram)
plotmat(t(transitMat))
pdf("transitionPlot.pdf")
plotmat(t(transitMat),box.size=0.05,shadow.size =0)
dev.off()

partTransitMat <- transitMat
partTransitMat[transitMat<0.1] <- NA
pdf("partTransitionPlot.pdf")
plotmat(t(partTransitMat),box.size=0.05,shadow.size =0)
dev.off()
}
