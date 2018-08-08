removeDeathOnFirstDay<-function(X,y){
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

X=X[lastDays>2,]
y=y[lastDays>2,]

n=nrow(X)
m=ncol(y)

X=cbind(rep(1,n),X)

ncovariates = ncol(X)
d = rep(-1,n)
for (j in 1:n) {
  deadDays = which(y[j,] == 0)
  if (length(deadDays) > 0)
    d[j] = min(deadDays)
}
lastDays <- d
lastDays[d == -1] = m
return(list(X=X,y=y,lastDays=lastDays))
}

initZ <- function(current_gamma)
  ###norders: total levels of y
{
  thresholds <- kappa(current_gamma)
  Z <- matrix(0,nr = n, nc = m, byrow = TRUE)
  indexZero <-
    y == 0   # we will set Z a random number below 0 for observations y==0.
  Z[indexZero] <-
    rtruncnorm(
      sum(indexZero), a = -Inf, b = 0, mean = -1, sd = 0.1
    )   
  indexSix <-
    y == (length(current_gamma) + 1)  # Set Z random numbers greater than kappa_5 for observations y==6.
  Z[indexSix] <-
    rtruncnorm(
      sum(indexSix), a = thresholds[length(current_gamma) + 2], b = Inf, mean = thresholds[length(current_gamma) +
                                                                                             2] + 1, sd = 0.1
    )
  for (i in 1:length(current_gamma))
  {
    # Set Z random numbers between kappa_{i-1} and kappa_i for observations y==i.
    indexMat <- y == i
    Z[indexMat] <-
      rtruncnorm(
        sum(indexMat), a = thresholds[i + 1], b = thresholds[i + 2], mean = 0.5 * (thresholds[i +
                                                                                                1] + thresholds[i + 2]), sd = 0.1
      )
  }
  Z
}

ul = function(y, current_gamma) {
  cutoffs<-kappa(current_gamma)
  return(c(cutoffs[y+2],cutoffs[y+1]))
}

log_pi_gamma = function(y, k,  hiMat, current_gamma, sigma2_epsilon) {
  n <- nrow(y)
  mu_gamma <- 0
  Sigma_gamma <- 10
  ret <- sum(log(dnorm(current_gamma,mu_gamma,Sigma_gamma)))  # log prior
  bound <- apply(y, c(1,2),ul,current_gamma)
  n_gamma <- length(current_gamma)
  indxMat <- (y>=k) & (y<=(n_gamma+1)) & (!is.na(hiMat))
  U_gamma <- bound[1,,]-hiMat
  L_gamma <- bound[2,,]-hiMat
  diff=pnorm(U_gamma[indxMat],sd=sqrt(sigma2_epsilon)) - pnorm(L_gamma[indxMat],sd=sqrt(sigma2_epsilon))
  return(ret+sum(log(diff[diff>0])))
}


logPI_gamma = function(y,  hiMat, current_gamma, sigma2_epsilon) {
  n <- nrow(y)
  mu_gamma <- 0
  Sigma_gamma <- 10
  ret <- sum(log(dnorm(current_gamma,mu_gamma,Sigma_gamma)))  # log prior
  bound <- apply(y, c(1,2),ul,current_gamma)
  n_gamma <- length(current_gamma)
  indxMat <-  (!is.na(hiMat))
  U_gamma <- bound[1,,]-hiMat
  L_gamma <- bound[2,,]-hiMat
  diff=pnorm(U_gamma[indxMat],sd=sqrt(sigma2_epsilon)) - pnorm(L_gamma[indxMat],sd=sqrt(sigma2_epsilon))
  return(ret+sum(log(diff[diff>0])))
}

kappa <- function(gamma)
{
  result<-rep(-1,length(gamma)+1)
  result[1]<-fixedThred1
  for(i in 2:(length(gamma)+1))
  {
    result[i]<-(result[i-1]+fixedThred2*exp(gamma[i-1]))/(1+exp(gamma[i-1]))
  }
 c(-Inf, result, fixedThred2, Inf)
}



zCategory <- function(current_gamma, zMat, y)
{
  thresholds <- kappa(current_gamma)
  zCat <- apply(zMat, c(1,2), findCat,thresholds)
  (zCat == y) * 1
}

findCat <- function(zvalue,thresholds)
{
  i <- 0
  while (thresholds[i + 1] < zvalue)
    i <- i + 1
  i-1
}

setKnots = function(yi, deadDay)
{
  lastDay = deadDay
  if (deadDay == -1)
    lastDay = m
  jumpInd=(yi[2:lastDay] -  yi[1:lastDay - 1]) != 0
  jumpPoints = sort(unique(c(1, (1:(lastDay-1))[jumpInd], lastDay)))
  return(jumpPoints)
}
