
oneStepGs <-
  function(bvec, sigma2_b,  cList, lambda, lambda_mu, gammaVec, zMat, beta, sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat) {
    accept <- rep(0,length(gammaVec))
    n <- nrow(X)
    a_b <- 1
    b_b <- 1/0.005
    scale_b <- (1 / b_b + 0.5 * sum((bvec)^2))
    post_a_b <- a_b  +0.5*n
    sigma2_b <- rinvgamma(1, shape = post_a_b, scale = scale_b)
    cRc_i_sum <- sum(sapply(1:n, function(i) t(cList[[i]]) %*% RmatList[[i]] %*% cList[[i]]))
    sumM <- sum(nbasisVec-2)
    b_lambda1 <- 1 / (1 / b_lambda + 0.5 * cRc_i_sum)    # scale
    lambda <-rgamma(1, shape = (sumM / 2 + a_lambda), scale = b_lambda1)*sigma2_epsilon
    for(i in 1:n)
    {
      sigma_c <-solve(basismat2List[[i]]  + lambda * RmatList[[i]])*sigma2_epsilon
      mu_c <-(sigma_c %*% t(basismatList[[i]]) %*% (zMat[i,1:lastObs[i]]-bvec[i]-matrix(X[i,]%*%beta,lastObs[i],1)))/sigma2_epsilon
      oneSamp <- mvrnorm(n = 1, mu_c, sigma_c)
      cList[[i]] <- oneSamp
    }
    muiMat <- t(sapply(1:n, function(i) c(basismatList[[i]] %*% cList[[i]], rep(NA, m-lastObs[i]))))
    mean_muiMat <- apply(muiMat,1,function(x) mean(x, na.rm=T))
    muiMat <- muiMat - matrix(rep(mean_muiMat,m), n,m,byrow = FALSE)
    post_sigma2_b <- 1/(m/sigma2_epsilon + 1/sigma2_b)
    post_mu_b <- post_sigma2_b*rowSums(zMat -muiMat- matrix(rep(X%*%beta,m),nr=n,nc=m),na.rm=TRUE)/sigma2_epsilon
    bvec <- rnorm(n, post_mu_b, sqrt(post_sigma2_b))
    for(i in 1:n){
      cList[[i]] <- cList[[i]]-mean_muiMat[i]
      if (lastObs[i] != m)
        cList[[i]][[nbasisVec[i]]] <- MIN - bvec[i]
    }
    rHiMat <- t(sapply(1:n, function(i) c(basismatList[[i]] %*% cList[[i]]+bvec[i], rep(MIN, m-lastObs[i]))))
    rHiMat[rHiMat>MAX] <- MAX
    rHiMat[rHiMat<MIN] <- MIN
    sigma_beta <- solve(sigma_0) + matrix(rowSums(apply(X, 1, function(x) {x%*%t(x)})%*%matrix(rep(m,n),nc=1)), nr=ncovariates)/sigma2_epsilon
    zMinusH <- zMat-rHiMat
    re <- apply(X,2,function(x) matrix(rep(x,m),nc=m)*zMinusH)
    mu_beta <- solve(sigma_0)%*%mu_0+colSums(re,na.rm = TRUE)/sigma2_epsilon
    sigma_beta <- solve(sigma_beta)
    mu_beta <- sigma_beta %*% mu_beta
    beta <- mvrnorm(n = 1, mu_beta, sigma_beta)
    hiMat <-  rHiMat+matrix(rep(X%*%beta,m),nr=n,nc=m)
    time <- proc.time()
    for (i in 1:length(gammaVec)) {
      log_pi_gamma_i <- function(gamma_i, i) {
        gammaVec1 <- gammaVec
        gammaVec1[i] <- gamma_i
        ret <- log_pi_gamma(y, i, hiMat, gammaVec1, sigma2_epsilon)
      }
      gamma_i_star <- rnorm(n = 1, gammaVec[i], sd_gamma_hat[i])
      log_mh <- log_pi_gamma_i(gamma_i_star, i) - log_pi_gamma_i(gammaVec[i],i)
      if (log_mh >= 0) {
        accept[i] <- 1
      }else{
        accept[i] <- exp(log_mh)
      }
      if (exp(log_mh) > runif(1,0,1)) {
        gammaVec[i] <- gamma_i_star
      }
    }
    proc.time()-time
    bound_gamma <- apply(y, c(1,2),ul,gammaVec)
    lowerBound <- bound_gamma[2,,]
    upperBound <- bound_gamma[1,,]
    zMat <- matrix(rtruncnorm(1, a = lowerBound, b = upperBound, mean = hiMat, sd = sqrt(sigma2_epsilon)),n,m)
    a_epsilon <- 1
    b_epsilon <- 1/0.005
    scale_epsilon <- (1 / b_epsilon + 0.5 * (sum((zMat -hiMat)^2,na.rm=T)))
    post_a_epsilon <- a_epsilon+0.5*n*m
    sigma2_epsilon <- rinvgamma(1, shape = post_a_epsilon, scale = scale_epsilon)
    list(bvec=bvec, sigma2_b=sigma2_b,  cList = cList,  lambda = lambda, lambda_mu = lambda_mu,   gammaVec = gammaVec, zMat = zMat, sigma2_epsilon =
           sigma2_epsilon, beta=beta, accept = accept
    )
  }

updateNiter <-
  function(niter, bvec, sigma2_b,  cList,  lambda, lambda_mu, gammaVec, zMat,beta,sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat)
  {
    oneGs <-
      oneStepGs(
        bvec, sigma2_b, cList, lambda, lambda_mu, gammaVec, zMat,beta, sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat
      )
    for (i in 1:niter) {
      cat("iter ", i,  "\n")
      oneGs <-
        oneStepGs(
          oneGs$bvec, oneGs$sigma2_b, oneGs$cList, oneGs$lambda, oneGs$lambda_mu, oneGs$gammaVec, oneGs$zMat, oneGs$beta, oneGs$sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat
        )
    }
    oneGs
  }

oneRunMCMC <-
  function(X,y, folder, nIterInOneBatch, nBatch, a_lambda, b_lambda, sd_gamma_hat) {
    if (!file.exists(folder)) {
      dir.create(folder)
    }   
    n <- nrow(X)
    filename_lambda = paste(folder,"lambda",sep = "")
    filename_beta = paste(folder,"beta",sep = "")
    filename_bvec = paste(folder,"bvec",sep = "")
    filename_sigma2_b = paste(folder,"sigma2_b",sep = "")
    filename_c = paste(folder,"c",sep = "")
    filename_gamma = paste(folder,"gamma",sep = "")
    filename_sigma2_epsilon = paste(folder,"sigma2_epsilon",sep = "")
    filename_accept = paste(folder,"accept",sep = "")
    method_gamma1 = 0
    method_gamma2 = 0
    cList = list()
    for (i in 1:n) {
      cList[[i]] = matrix(rep(0.5,nbasisVec[i]),ncol=1)
    }
    for(i in 1:n) if(lastDays[i]!=m) cList[[i]][nbasisVec[i]]=MIN
    lambda_mu <- 35
    lambda <- 35
    gammaVec <- rnorm(4)
    zMat <-  matrix(rnorm(n*m),n,m) # initZ(gammaVec)
    beta <- c(4, 0.2,-0.08)+rep(0, ncol(X))
    sigma2_epsilon=0.5
    sigma2_b <- 1
    bvec <-  rnorm(n,0, sqrt(sigma2_b))
    ifAppend = FALSE;
      gsRe <-
        oneStepGs(
          bvec, sigma2_b, cList, lambda, lambda_mu, gammaVec, zMat,beta, sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat
        )
    for (j in 1:nBatch) {
      cat("batach: ", j, "\n")
        
        for (i in 1:nIterInOneBatch) {
          cat("iter ", i,  "\n")
          gsRe <-
            oneStepGs(
              gsRe$bvec, gsRe$sigma2_b, gsRe$cList, gsRe$lambda, gsRe$lambda_mu, gsRe$gammaVec, gsRe$zMat, gsRe$beta, gsRe$sigma2_epsilon, a_lambda, b_lambda, sd_gamma_hat
            )
        }
      
      if(j>2)ifAppend = TRUE;
      write.table(
        #    matrix(gsRe$lambdaMat,nr = 1), file = paste(filename_lambda,".txt",sep =
        gsRe$lambda, file = paste(filename_lambda,".txt",sep =
                                    ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );
      
      write.table(
        matrix(gsRe$bvec,nrow=1), file = paste(filename_bvec,".txt",sep =
                                                 ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );
      
      write.table(
        gsRe$sigma2_b, file = paste(filename_sigma2_b,".txt",sep =
                                      ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );
      
      write.table(
        matrix(gsRe$gamma,nrow=1), file = paste(filename_gamma,".txt",sep = ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );
      
      write.table(
        matrix(gsRe$beta,nrow=1), file = paste(filename_beta,".txt",sep = ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );

      write.table(
        matrix(gsRe$sigma2_epsilon,ncol=1), file = paste(filename_sigma2_epsilon,".txt",sep = ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );
      
      write.table(
        matrix(gsRe$accept,nrow=1), file = paste(filename_accept,".txt",sep = ""), append = ifAppend, row.names = FALSE,
        col.names = FALSE
      );

      for (i in 1:n)
      {
        cvec <- matrix(gsRe$cList[[i]],nrow=1)
        write.table(
          cvec, file = paste(filename_c,"_",i, ".txt", sep = ""), append = ifAppend, row.names = FALSE, col.names = FALSE
        );
      }
    }
    return(basismatList)
  }

