##--------------------------------------------------##
##           One-shot simulation study              ##
##           of multivariate regression             ##
##        with element-wise contamination           ##
##--------------------------------------------------##

## Preliminary
rm(list=ls())
library(MASS)

source("CSM-Reg.R")
source("MVT-Reg.R")

## settings
p <- 6     # dimension of y
qq <- 10   # number of covariates 
n <- 200
contami_ratio <- 0.5

lab <- c("CSM", "CT", "AT", "G")


## true parameters
Beta_true <- rep(0, qq)
Beta_true[1] <- 0.5 
Beta_true[2] <- 1
Beta_true[3] <- (-1) 

Sig_true <- diag(p)  
for(i in 1:p){
  for(j in 1:p){
    Sig_true[i,j] <- (0.6)^(abs(i-j)) 
  }
}


## data generation 
mat <- 0.7*diag(qq) + 0.3*matrix(1, qq, qq)
X <- array(NA, c(n, p, qq))
for(i in 1:n){
  X[i,,] <- mvrnorm(p, rep(0,qq), mat)
}

Y <- matrix(NA, n, p)
for(i in 1:n){
  Mu <- as.vector( X[i,,]%*%Beta_true )
  Y[i,] <- mvrnorm(1, Mu, Sig_true)
}


## outlier
Z <- matrix(0, n, p)
for(i in 1:n){
  if( runif(1)<contami_ratio ){
    sel <- sample(1:p, 1)
    Y[i, sel] <-  Y[i, sel] + 10
    Z[i, sel] <- 1
  }
}


## proposed CSM method
mc <- 2000
bn <- 1000
fit_CSM <- CSM_Reg(Y, X, mc=mc, bn=bn)

## t-distribution 
fit_CT <- CT_Reg(Y, X, mc=mc, bn=bn)
fit_AT <- AT_Reg(Y, X, mc=mc, bn=bn)

## Gaussian 
fit_G <- CT_Reg(Y, X, mc=mc, bn=bn, t_dist=F)


## Regression coefficients
hBeta_CSM <- apply(fit_CSM$Beta, 2, mean)
hBeta_CT <- apply(fit_CT$Beta, 2, mean)
hBeta_AT <- apply(fit_AT$Beta, 2, mean)
hBeta_G <- apply(fit_G$Beta, 2, mean)

Est <- cbind(hBeta_CSM, hBeta_CT, hBeta_AT, hBeta_G)
dimnames(Est)[[2]] <- lab
100*sqrt( apply((Est-Beta_true)^2, 2, mean) )   # estimation error


## Error covariance matrix
hSig_CSM <- apply(fit_CSM$Sig, c(2,3), mean)
hSig_CT <- apply(fit_CT$Sig, c(2,3), mean)
hSig_AT <- apply(fit_AT$Sig, c(2,3), mean)
hSig_G <- apply(fit_G$Sig, c(2,3), mean)

MSE_sig <- rep(NA, 4)
names(MSE_sig) <- lab
MSE_sig[1] <- mean((hSig_CSM-Sig_true)^2)    
MSE_sig[2] <- mean((hSig_CT-Sig_true)^2)
MSE_sig[3] <- mean((hSig_AT-Sig_true)^2)
MSE_sig[4] <- mean((hSig_G-Sig_true)^2)
MSE_sig 


## Credible intervals of regression coefficients 
CI <- list()
CI[[1]] <- apply(fit_CSM$Beta, 2, quantile, prob=c(0.025, 0.975))
CI[[2]] <- apply(fit_CT$Beta, 2, quantile, prob=c(0.025, 0.975))
CI[[3]] <- apply(fit_AT$Beta, 2, quantile, prob=c(0.025, 0.975))
CI[[4]] <- apply(fit_G$Beta, 2, quantile, prob=c(0.025, 0.975))


CP <- AL <- rep(NA, 4)
names(CP) <- names(AL) <- c("CSM", "CT", "AT", "G")
for(l in 1:4){
  CP[l] <- mean(CI[[l]][1,]<Beta_true & CI[[l]][2,]>Beta_true)    # coverage probability
  AL[l] <- mean(CI[[l]][2,]-CI[[l]][1,])      # average interval length
}

CP
AL
