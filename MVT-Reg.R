library(MASS)
library(MCMCpack)
library(progress)


###-------------------------------------------###
###       Multivariate regression with        ###
###       classical t-distribution (CT)       ###
###-------------------------------------------###
# Y: n*p-matrix

CT_Reg <- function(Y, X, nu=5, nu_est=T, mc=2000, bn=1000, t_dist=T){
  # preliminary 
  p <- ncol(Y)
  n <- nrow(Y)
  qq <- dim(X)[3]
  n0 <- p
  S0 <- (n0+p+1)*diag(p)
  mu <- rep(0, qq)
  Psi <- 0.01*diag(qq)
  Psi_mu <- c(Psi%*%mu)
  
  # Initial values
  Beta <- rep(0, qq)    # Initialize regression coefficient
  Sigma <- var(Y)  # Initialize covariance matrix
  lam <- rep(1, n)  # Initialize scale parameter
  
  # Matrices to store posterior samples
  Beta_pos <- matrix(NA, mc, qq)
  Sigma_pos <- array(0, c(mc, p, p))
  nu_pos <- rep(NA, mc)
  
  # Gibbs sampling
  pb <- progress_bar$new(total=mc)   # progress bar 
  for(iter in 1:mc){
    # residual 
    E <- matrix(NA, n, p)
    for(i in 1:n){
      E[i,] <- Y[i,] - c(X[i,,]%*%Beta)
    }
    
    # lambda
    if(t_dist){
      for(i in 1:n){
        shape <- (nu + p) / 2
        scale <- (nu + t(E[i,])%*%solve(Sigma)%*%E[i,]) / 2
        lam[i] <- 1 / rgamma(1, shape, scale)
      }
    }
    
    # nu (degrees of freedom)
    if(t_dist & nu_est){
      nu_new <- nu + rnorm(1)
      if(nu_new<1){ nu_new <- 1 }
      if(nu_new>100){ nu_new <- 100 }
      log_dens <- sum( dgamma(1/lam, 0.5*nu, 0.5*nu, log=T) )
      log_dens_new <- sum( dgamma(1/lam, 0.5*nu_new, 0.5*nu_new, log=T) )
      pp <- exp(log_dens_new - log_dens)
      if(runif(1)<pp){
        nu <- nu_new
      }
    }
    
    # Sigma
    S <- matrix(0, p, p)
    for (i in 1:n) {
      S <- S + E[i,]%*%t(E[i,])/lam[i]
    }
    Sigma <- riwish(n0+n, S0+S)
    inv_Sigma <- solve(Sigma)
    
    # Beta
    At <- sapply(1:n, function(i){
      X_i <- X[i,,]
      return( t(X_i)%*%inv_Sigma%*%X_i/lam[i] )
    })
    A <- matrix(rowSums(At), qq, qq)
    bec <- sapply(1:n, function(i){
      X_i <- X[i,,]
      return( c(t(X_i)%*%(inv_Sigma%*%(Y[i,]/lam[i]))) )
    })
    bec <- rowSums(bec)
    solve_Precision <- solve(A + Psi)
    Beta <- mvrnorm(1, mu=c(solve_Precision%*%(Psi_mu+bec)), Sigma=solve_Precision)
    
    # Store samples
    Beta_pos[iter,] <- Beta
    Sigma_pos[iter,,] <- Sigma
    nu_pos[iter] <- nu
    
    # progress
    pb$tick()
  }
  
  # Excluding burn-in period
  Sigma_pos <- Sigma_pos[-(1:bn),,]
  Beta_pos <- Beta_pos[-(1:bn),]
  nu_pos <- nu_pos[-(1:bn)]
  
  # return
  Result <- list(Beta=Beta_pos, Sigma=Sigma_pos, nu=nu_pos)
  return(Result)
}








###-------------------------------------------###
###       Multivariate regression with        ###
###      alternative t-distribution (AT)      ###
###-------------------------------------------###
# Y: n*p-matrix

AT_Reg <- function(Y, X, nu=5, nu_est=T, mc=2000, bn=1000){
  # preliminary 
  p <- ncol(Y)
  n <- nrow(Y)
  qq <- dim(X)[3]
  n0 <- p
  S0 <- (n0+p+1)*diag(p)
  mu <- rep(0, qq)
  Psi <- 0.01*diag(qq)
  Psi_mu <- c(Psi%*%mu)
  
  # Initial values
  Beta <- rep(0, qq)    # Initialize regression coefficient
  Sigma <- var(Y)  # Initialize covariance matrix
  tau <- matrix(1, n, p)  # Initialize scale parameter
  
  # Matrices to store posterior samples
  Beta_pos <- matrix(NA, mc, qq)
  Sigma_pos <- array(0, c(mc, p, p))
  nu_pos <- rep(NA, mc)
  
  # Gibbs sampling
  pb <- progress_bar$new(total=mc)   # progress bar 
  for(iter in 1:mc){
    # residual 
    E <- matrix(NA, n, p)
    for(i in 1:n){
      E[i,] <- Y[i,] - c(X[i,,]%*%Beta)
    }
    
    # tau
    Om <- solve(Sigma)
    for(i in 1:n){
      for(k in 1:p){
        tau_new <- tau[i,k] + 0.2*rnorm(1)
        if(tau_new<10^(-3)){  tau_new <- 10^(-3)  }
        alpha <- 0.5*(nu+1)
        be <- 0.5*nu + 0.5*E[i,k]^2*Om[k,k]
        U_ik <- E[i,-k]*sqrt(tau[i,-k])
        gam <- E[i,k]*sum(Om[-k,k]*U_ik)
        log_dens <- (alpha-1)*log(tau[i,k]) - be*tau[i,k] - gam*sqrt(tau[i,k])
        log_dens_new <- (alpha-1)*log(tau_new) - be*tau_new - gam*sqrt(tau_new)
        pp <- exp(log_dens_new - log_dens)
        if(runif(1)<pp){
          tau[i,k] <- tau_new
        }
      }
    }
    
    # nu (degrees of freedom)
    if(nu_est){
      nu_new <- nu + rnorm(1)
      if(nu_new<1){ nu_new <- 1 }
      if(nu_new>100){ nu_new <- 100 }
      log_dens <- sum( dgamma(tau, 0.5*nu, 0.5*nu, log=T) )
      log_dens_new <- sum( dgamma(tau, 0.5*nu_new, 0.5*nu_new, log=T) )
      pp <- exp(log_dens_new - log_dens)
      if(runif(1)<pp){
        nu <- nu_new
      }
    }
    
    # Sigma
    S <- matrix(0, p, p)
    for(i in 1:n){
      E_vec <- E[i,]*sqrt(tau[i,])
      S <- S + E_vec%*%t(E_vec)
    }
    Sigma <- riwish(n0+n, S0+S)
    inv_Sigma <- solve(Sigma)
    
    # Beta
    At <- sapply(1:n, function(i){
      X_i <- X[i,,]
      Xt_i <- sqrt(tau[i,])*X_i
      return( t(Xt_i)%*%inv_Sigma%*%Xt_i )
    })
    A <- matrix(rowSums(At), qq, qq)
    bec <- sapply(1:n, function(i){
      X_i <- X[i,,]
      Xt_i <- sqrt(tau[i,])*X_i
      return( c(t(Xt_i)%*%(inv_Sigma%*%(Y[i, ]*sqrt(tau[i,])))) )
    })
    bec <- rowSums(bec)
    solve_Precision <- solve(A + Psi)
    Beta <- mvrnorm(1, mu=c(solve_Precision%*%(Psi_mu+bec)), Sigma=solve_Precision)
    
    # Store samples
    Beta_pos[iter,] <- Beta
    Sigma_pos[iter,,] <- Sigma
    nu_pos[iter] <- nu
    
    # progress
    pb$tick()
  }
  
  # Excluding burn-in period
  Sigma_pos <- Sigma_pos[-(1:bn),,]
  Beta_pos <- Beta_pos[-(1:bn),]
  nu_pos <- nu_pos[-(1:bn)]
  
  # return
  Result <- list(Beta=Beta_pos, Sigma=Sigma_pos, nu=nu_pos)
  return(Result)
}

