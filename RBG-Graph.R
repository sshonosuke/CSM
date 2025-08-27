# Bayesian Robust sparce covariance estimation

# library for glasso
library(glasso)
library(glassoFast)
library(MCMCpack)



RBGGM.function <- function(Y, gam=0.1, mc=1000, lambda=1, core=1){
  
  library(parallel)
  DC <- detectCores()
  if(DC<= core){
    core <- DC
  }
  
  fit0 <- FRGGM.function(Y, gam, eps=10^(-5), lambda)
  
  set.seed(100)
  
  fun  <- function(i){
    #set.seed(i)
    
    p <- dim(Y)[2]
    n <- dim(Y)[1]
    
    eps <- 10^(-3)
    
    Ome.pos <- matrix(NA, p,p)
    Fgam <- rep(NA,n)
    
    Sig <- t(Y)%*%(Y)/n
    
    expbox <- rep(NA,n)
    ## Generate weight from Dirichlet distribution
    ww <- (n+1)*as.vector(rdirichlet(1,rep(1,n+1)))
    w <- ww[1:n]
    wp <- ww[n+1]
    if(wp*lambda < 0.001){wp <- 0.001/lambda}
    
    ## Sample approximate posterior samples via WBB
    LAM <- 2*(1+gam)*wp*lambda
    
    oldOme <- Ome <- Omef <- fit0
    Sigma <- solve(fit0)
    
    for (k in 1:10000) {
      
      ## setting 
      for (j in 1:n) {
        expbox[j] <- exp(-gam*t(Y[j,])%*%Ome%*%Y[j,]/2)
        Fgam[j] <- (2*pi)^(-p*gam/2)*det(Ome)^(gam/2)*expbox[j]
      }
      wF <- w*Fgam
      sumwF <- sum(wF)
      s <- as.vector(wF/sumwF)
      
      S <- 0
      for (j in 1:n) {
        S <- S + s[j]*(Y[j,])%*%t(Y[j,])
      }
      S <- (1+gam)*(S+t(S))/2
      
      oldOme <- Ome
      
      ## graphical Lasso
      fit <- glassoFast(S=S, rho=LAM, wi.init= oldOme)
      Ome <- fit$wi
 
      if(max(abs(Ome-oldOme)) < eps){
        break
      }
    }
    
    return(Ome)
  }
  
  d <- mclapply(1:mc, fun, mc.cores=core)
  DFfit <- array(0, dim=c(p,p,mc))
  for(i in 1:mc){
    DFfit[,,i] <- d[[i]]
  }
  DF <- DFfit
  
  print(paste0("WBB samples: ",mc))
  return(DF)
}







FRGGM.function <- function(Y, gam, eps=10^(-5), lambda=1){
  
  p <- dim(Y)[2]
  n <- dim(Y)[1]
  
  Sig <- t(Y)%*%(Y)/n
  Fgam <- rep(NA,n)
  LAM <- 2*(1+gam)*lambda
  
  fit <- glassoFast(S=Sig, rho=LAM)
  Ome <- fit$wi
  
  for (m in 1:10000) {

    for (i in 1:n) {
      Fgam[i] <- ( (2*pi)^(-p/2)*det(Ome)^(1/2)*exp(-t(Y[i,])%*%Ome%*%Y[i,]/2) )^(gam)
    }
    sumF <- sum(Fgam)
    s <- as.vector(Fgam/sumF)
    S <- 0
    for (j in 1:n) {
      S <- S + s[j]*(Y[j,])%*%t(Y[j,])
    }
    S <- (1+gam)*(S+t(S))/2
    
    oldOme <- Ome
    fit <- glassoFast(S=S, rho=LAM, wi.init=oldOme)
    Ome <- fit$wi
    
    if(max(abs(Ome-oldOme)) < eps){
      break
    }
  }

  return(Ome)
}
