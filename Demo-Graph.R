##--------------------------------------------------##
##           One-shot simulation study              ##
##            of graphical modeling                 ##
##        with element-wise contamination           ##
##--------------------------------------------------##

## Preliminary
rm(list=ls())
library(MASS)
source("CSM-Graph.R")   # proposed CSM method 
source("RBG-Graph.R")   # gamma divergence 
source("MVT-Graph.R")   # t-distribution 


## settings
scenario <- 1           #  1-3
contami_ratio <- 0.2    #  0.2, 0.4, 0.6
p <- 12
n <- 200
set.seed(scenario)


## true precision matrix 
Omega <- diag(p)  
for(i in 1:p){
  for(j in 1:p){
    if(abs(i-j)==1){ Omega[i,j] <- 0.5 }
    if(abs(i-j)==2){ Omega[i,j] <- 0.25 }
  }
}


## data generation 
Y <- mvrnorm(n, rep(0,p), solve(Omega))

Z <- matrix(0, n, p)  # outlier indicator 
for(i in 1:n){
  if( runif(1)<contami_ratio ){
    if(scenario<3){
      sel <- sample(1:p, scenario)
      Y[i, sel] <-  Y[i, sel] + 10
      Z[i, sel] <- 1
    }
    if(scenario==3){
      num <- 1+rpois(1, lambda=1)
      sel <- sample(1:p, num)
      Y[i, sel] <-  Y[i, sel] + 10
      Z[i, sel] <- 1
    }
  }
}


## proposed CSM method
mc <- 2000
bn <- 1000
fit_CSM <- CSM_Graph(Y, mc=mc, bn=bn) 
hOm_CSM <- apply(fit_CSM$Om, c(1,2), mean)

## Bayesian gamma-divergence
fit_RGB1 <- RBGGM.function(Y, gam=0.05, mc=1000, lambda=0.01)
hOm_RGB1 <- apply(fit_RGB1, c(1,2), mean)
fit_RGB2 <- RBGGM.function(Y, gam=0.1, mc=1000, lambda=0.01)
hOm_RGB2 <- apply(fit_RGB2, c(1,2), mean)

## t-distribution
fit_cT <- MVT(Y)
hOm_cT <- apply(fit_cT$Om, c(2,3), mean)
fit_aT <- AMVT(Y)
hOm_aT <- apply(fit_aT$Om, c(2,3), mean)
fit_dT <- DMVT(Y)
hOm_dT <- apply(fit_dT$Om, c(2,3), mean)

## Gaussian
fit_GG <- MVT(Y, mc=mc, bn=bn, t_dist=F)
hOm_GG <- apply(fit_GG$Om, c(2,3), mean)



## estimation error 
mean((Omega-hOm_CSM)^2)    # proposed CSM method 
mean((Omega-hOm_RGB1)^2)   # Bayesian gamma divergence (gamma=0.05)
mean((Omega-hOm_RGB2)^2)   # Bayesian gamma divergence (gamma=0.1)
mean((Omega-hOm_cT)^2)     # classical t-distribution 
mean((Omega-hOm_aT)^2)     # alternative t-distribution 
mean((Omega-hOm_dT)^2)     # Dirichlet t-distribution 
mean((Omega-hOm_GG)^2)     # Gaussian 

