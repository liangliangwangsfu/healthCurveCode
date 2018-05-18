m = ncol(y)
M = m - 2
n = nrow(y)  # number of patients

setKnots = function(yi, deadDay)
{
  lastDay = deadDay
  if (deadDay == -1)
    lastDay = m
  jumpInd=(yi[2:lastDay] -  yi[1:lastDay - 1]) != 0
  jumpPoints = sort(unique(c(1, (1:(lastDay-1))[jumpInd], lastDay)))
  return(jumpPoints)
}

#basisFunctions<-function(m,n){
# number of quadrature points (including two knots) between two knots
norder   = 4
tobsList = list()
knotsList = list()
nbasisVec = rep(0,n)
bsbasisList=list()
basismatList = list()
basismat2List = list()
DbasismatList = list()
D2basismatList = list()
RmatList = list()
basismatAllDaysList=list()
lastObs=lastDays
if(knotsOption==1)
{
  K = m+2
  commonKnots=seq(1,m,length.out = K-2)  
}  else
  if(knotsOption==2)
  {
    commonKnots=seq(1,m,by=distBetween2knots)
  } 
  
hatList = list()
for (i in 1:n) {
  tobsList[[i]] = 1:lastObs[i]
  if(knotsOption==3)
  knotsList[[i]] = setKnots(y[i,], lastObs[i]) #c(commonKnots[commonKnots<lastObs[i]],lastObs[i])
  else
  knotsList[[i]] = c(commonKnots[commonKnots<lastObs[i]],lastObs[i])
  nbasisVec[i]   =   length(knotsList[[i]])  + norder - 2
  bsbasisList[[i]]  = create.bspline.basis(range(knotsList[[i]]),nbasisVec[i],norder,knotsList[[i]])
  # basis values at sampling points
  basismatList[[i]]   = eval.basis(tobsList[[i]],   bsbasisList[[i]])
  # square of basis matrix (Phi). It is used frequent for optimization,
  # so we calculate in advance
  basismat2List[[i]]  = t(basismatList[[i]]) %*% basismatList[[i]]
  # values of the first derivative of basis functions at sampling points
  DbasismatList[[i]]  = eval.basis(tobsList[[i]],  bsbasisList[[i]], 1)
  # values of the second derivative of basis functions at sampling points
  D2basismatList[[i]] = eval.basis(tobsList[[i]],   bsbasisList[[i]], 2)
  nquad                        = 5;
  # set up quadrature points and weights
  #[sinbasis, quadpts, quadwts]
  quadre = quadset(nquad, bsbasisList[[i]])
  quadpts = quadre$quadvals[,1]
  quadwts = quadre$quadvals[,2]
  # values of the second derivative of basis functions at quadrature points
  D2quadbasismat = eval.basis(quadpts, bsbasisList[[i]], 2)
  RmatList[[i]]  = t(D2quadbasismat) %*% (D2quadbasismat * ((quadwts) %*%
                                                              matrix(1,1,nbasisVec[i])))
}

