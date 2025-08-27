library(MASS)
library(MCMCpack)
library(progress)


###-------------------------------------------###
###         Graphical modeling with           ###
###         classical t-distribution          ###
###-------------------------------------------###
# Y: n*p-matrix

MVT <- function(Y, nu=5, nu_est=T, mc=2000, bn=1000, t_dist=T){
  # preliminary 
  n <- nrow(Y)
  p <- ncol(Y)
  n0 <- p
  S0 <- (2*p+1)*diag(p)
  
  # Initial values
  Sigma <- var(Y)  # Initialize covariance matrix
  lam <- rep(1, n)  # Initialize scale parameter
  
  # Matrices to store posterior samples
  Sigma_pos <- array(0, c(mc, p, p))
  Om_pos <- array(0, c(mc, p, p))
  
  # Gibbs sampling
  pb <- progress_bar$new(total=mc)   # progress bar 
  for(iter in 1:mc){
    # lambda
    if(t_dist){
      for(i in 1:n){
        shape <- (nu + p) / 2
        scale <- (nu + t(Y[i,])%*%solve(Sigma)%*%(Y[i,])) / 2
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
      S <- S + Y[i,]%*%t(Y[i,])/lam[i]
    }
    Sigma <- riwish(n0+n, S0+S)
    
    # Store samples
    Sigma_pos[iter,,] <- Sigma
    Om_pos[iter,,] <- solve(Sigma)
    
    pb$tick()
  }
  
  # Check results (excluding burn-in period)
  Sigma_pos <- Sigma_pos[-(1:bn),,]
  Om_pos <- Om_pos[-(1:bn),,]
  
  # return
  Result <- list(Sigma=Sigma_pos, Om=Om_pos)
  return(Result)
}





###-------------------------------------------###
###         Graphical modeling with           ###
###       alternative t-distribution          ###
###-------------------------------------------###
# Y: n*p-matrix
AMVT <- function(Y, nu=5, nu_est=T, mc=2000, bn=1000){
  # preliminary 
  n <- nrow(Y)
  p <- ncol(Y)
  n0 <- p
  S0 <- (2*p+1)*diag(p)
  
  # Initial values
  mu <- rep(0, p)
  Sigma <- var(Y)  # Initialize covariance matrix
  tau <- matrix(1, n, p)  # Initialize scale parameter
  
  # Matrices to store posterior samples
  Sigma_pos <- array(0, c(mc, p, p))
  Om_pos <- array(0, c(mc, p, p))
  nu_pos <- rep(0, mc)
  Tau_pos <- array(0, c(mc, n, p))
  
  # Gibbs sampling
  pb <- progress_bar$new(total=mc)   # progress bar 
  for(iter in 1:mc){
    # tau
    Om <- solve(Sigma)
    for(i in 1:n){
      for(k in 1:p){
        tau_new <- tau[i,k] + 0.2*rnorm(1)
        if(tau_new<10^(-3)){  tau_new <- 10^(-3)  }
        alpha <- 0.5*(nu+1)
        be <- 0.5*nu + 0.5*Y[i,k]^2*Om[k,k]
        X_ik <- Y[i,-k]*sqrt(tau[i,-k])
        gam <- Y[i,k]*sum(Om[-k,k]*X_ik)
        log_dens <- (alpha-1)*log(tau[i,k]) - be*tau[i,k] - gam*sqrt(tau[i,k])
        log_dens_new <- (alpha-1)*log(tau_new) - be*tau_new - gam*sqrt(tau_new)
        pp <- exp(log_dens_new - log_dens)
        if(runif(1)<pp){
          tau[i,k] <- tau_new
        }
      }
    }
    Tau_pos[iter,,] <- tau
    
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
      X_vec <- Y[i,]*sqrt(tau[i,])
      S <- S + X_vec%*%t(X_vec)
    }
    Sigma <- riwish(n0+n, S0+S)
    
    # Store samples
    Sigma_pos[iter,,] <- Sigma
    Om_pos[iter,,] <- solve(Sigma)
    nu_pos[iter] <- nu
    
    pb$tick()
  }
  
  # Check results (excluding burn-in period)
  Sigma_pos <- Sigma_pos[-(1:bn),,]
  Om_pos <- Om_pos[-(1:bn),,]
  nu_pos <- nu_pos[-(1:bn)]
  Tau_pos <- Tau_pos[-(1:bn),,]
  
  # return
  Result <- list(Sigma=Sigma_pos, Om=Om_pos, nu=nu_pos, Tau=Tau_pos)
  return(Result)
}







###-------------------------------------------###
###         Graphical modeling with           ###
###         Dirichlet t-distribution          ###
###-------------------------------------------###
# Y: n*p-matrix

DMVT <- function(Y, nu=3, nu_est=F, mc=2000, bn=1000){
  # preliminary 
  n <- nrow(Y)
  p <- ncol(Y)
  n0 <- p
  S0 <- (2*p+1)*diag(p)
  
  # Initial values
  Sigma <- var(Y)  # Initialize covariance matrix
  K <- round(0.8*p)
  tau_base <- matrix(1, n, K)  # (group-wise) scale parameter
  for(i in 1:n){
    tau_base[i,] <- seq(1, 0, length=K+1)[1:K]
  }
  Z <- matrix(1, n, p)    # group assignment 
  Tau <- matrix(1, n, p)   # scale parameter
  Pi <- rep(1/K, K)
  d_alpha <- 3       # parameter of Dirichlet distribution 
  
  # Matrices to store posterior samples
  Sigma_pos <- array(0, c(mc, p, p))
  Om_pos <- array(0, c(mc, p, p))
  nu_pos <- rep(0, mc)
  Tau_pos <- array(0, c(mc, n, p))
  
  # Gibbs sampling
  pb <- progress_bar$new(total=mc)   # progress bar 
  for(iter in 1:mc){
    # Z
    Om <- solve(Sigma)
    for(i in 1:n){
      for(j in 1:p){
        alpha <- 0.5
        be <- 0.5*Y[i,j]^2*Om[j,j]
        X_ij <- Y[i,-j]*sqrt(Tau[i,-j])
        gam <- Y[i,j]*sum(Om[-j,j]*X_ij)
        log_dens <- log(Pi+10^(-5)) + (alpha-1)*log(tau_base[i,]) - be*tau_base[i,] - gam*sqrt(tau_base[i,])
        log_dens <- log_dens - max(log_dens)
        prob <- exp(log_dens) / sum(exp(log_dens))
        Z[i,j] <- sample(1:K, 1, prob=prob)
      }
    }
    
    # Tau (base)
    for(i in 1:n){
      for(k in 1:K){
        Nk <- sum(Z[i,]==k)
        if(Nk>0){
          tau_new <- tau_base[i,k] + 0.2*rnorm(1)
          if(tau_new<10^(-3)){  tau_new <- 10^(-3)  }
          alpha <- 0.5*nu + 0.5*Nk
          sY <- Y[i,Z[i,]==k]
          rY <- Y[i,Z[i,]!=k]
          sOm <- diag( Om[Z[i,]==k, Z[i,]==k] )
          rOm <- Om[Z[i,]==k, Z[i,]!=k] 
          be <- 0.5*nu + 0.5*sum(sY^2*sOm)
          X_ij <- rY*sqrt(Tau[i,Z[i,]!=k])
          if(length(X_ij)>1){
            gam <- sum( sY*c(rOm%*%X_ij) )
          }else{
            gam <- sum( sY*rOm*X_ij )
          }
          log_dens <- (alpha-1)*log(tau_base[i,k]) - be*tau_base[i,k] - gam*sqrt(tau_base[i,k])
          log_dens_new <- (alpha-1)*log(tau_new) - be*tau_new - gam*sqrt(tau_new)
          pp <- exp(log_dens_new - log_dens)
          if(runif(1)<pp){
            tau_base[i,k] <- tau_new
          }
        }else{
          tau_base[i,k] <- rgamma(1, 0.5*nu, 0.5*nu)
        }
        Tau[i, Z[i,]==k] <- tau_base[i,k]
      }
    }
    Tau_pos[iter,,] <- Tau
    
    # re-labeling
    new_Z <- Z
    new_tau_base <- tau_base
    for(i in 1:n){
      ord <- order(as.numeric(table( factor(Z[i,], levels=1:K) )), decreasing=T)
      for(k in 1:K){
        new_Z[i, Z[i,]==ord[k]] <- k
        new_tau_base[i, k] <- tau_base[i, ord[k]]
      }
    }
    Z <- new_Z
    tau_base <- new_tau_base
    
    # Pi
    alpha_t <- rep(d_alpha, K) + as.numeric(table( factor(Z, levels=1:K) ))
    Pi <- c(rdirichlet(1, alpha_t))
    
    # d_alpha
    d_alpha_new <- d_alpha + 0.3*rnorm(1)
    if(d_alpha_new<0.1){ d_alpha_new <- 0.1 }
    if(d_alpha_new>100){ d_alpha_new <- 100 }
    log_dens <- dgamma(d_alpha, 1, 1, log=T) + log(ddirichlet(Pi, rep(d_alpha, K))) 
    log_dens_new <- dgamma(d_alpha_new, 1, 1, log=T) + log(ddirichlet(Pi, rep(d_alpha_new, K))) 
    pp <- exp(log_dens_new - log_dens)
    if(runif(1)<pp){
      d_alpha <- d_alpha_new
    }
    
    # nu (degrees of freedom)
    if(nu_est){
      nu_new <- nu + rnorm(1)
      if(nu_new<1){ nu_new <- 1 }
      if(nu_new>100){ nu_new <- 100 }
      log_dens <- sum( dgamma(tau_base, 0.5*nu, 0.5*nu, log=T) )
      log_dens_new <- sum( dgamma(tau_base, 0.5*nu_new, 0.5*nu_new, log=T) )
      pp <- exp(log_dens_new - log_dens)
      if(runif(1)<pp){
        nu <- nu_new
      }
    }
    
    # Sigma
    S <- matrix(0, p, p)
    for(i in 1:n){
      X_vec <- Y[i,]*sqrt(Tau[i,])
      S <- S + X_vec%*%t(X_vec)
    }
    Sigma <- riwish(n0+n, S0+S)
    
    # Store samples
    Sigma_pos[iter,,] <- Sigma
    Om_pos[iter,,] <- solve(Sigma)
    nu_pos[iter] <- nu
    
    pb$tick()
  }
  
  # Check results (excluding burn-in period)
  Sigma_pos <- Sigma_pos[-(1:bn),,]
  Om_pos <- Om_pos[-(1:bn),,]
  nu_pos <- nu_pos[-(1:bn)]
  Tau_pos <- Tau_pos[-(1:bn),,]
  
  # return
  Result <- list(Sigma=Sigma_pos, Om=Om_pos, nu=nu_pos, Tau=Tau_pos)
  return(Result)
}









