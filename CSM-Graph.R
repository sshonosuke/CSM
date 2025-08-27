library(MASS)
library(progress)
#install.packages("tmg_0.3.tar", repos=NULL, type="source")  # for installing "tmg" package
library(tmg)




###----------------------------------------------------###
###        Robust Gaussian Graphical Modeling          ###
###     using Correlation-intact Sandwich Mixture      ###
###----------------------------------------------------###
# Y: n*p-matrix  (Y should be centered) 

CSM_Graph <- function(Y, mc=2000, bn=1000, positive=F){
  # preliminary
  upper_bound <- 10^8
  p <- ncol(Y)
  n <- nrow(Y)
  
  # hyper-parameters of priors 
  de <- 10^(-8)
  de_dash <- 0
  de_tilde <- 1
  ga <- 1
  a <- 0.05 
  b <- 1  
  nu <- p
  Ga <- (nu+p+1)*diag(p)
  
  # function (rejection)
  rejection <- function(ga, de){
    TF <- FALSE
    while(!TF){
      r <- rbeta(1, shape1=1, shape2=ga)
      r <- exp( r/(1-r) )
      TF <- ( log(runif(1)) <= (-de)/(2*r^2) )
    }
    return(r)
  }
  
  # initial values 
  r <- matrix(2, n, p)    # latent variable
  s <- matrix(1, n, p)    # latent variable
  z <- matrix(1, n, p)    # latent variable
  t. <- r*s
  tz <- t.^z
  phi <- 0.05    # outlier probability
  Si <- diag(p)
  solve_Si <- solve(Si)
  psi <- sqrt( diag(solve_Si) )
  Q <- t( t(solve_Si)/psi )
  Q <- Q/psi
  
  # objects for posterior samples
  Si_pos <- array(NA, c(p, p, mc))
  z_pos <- array(NA, c(n, p, mc))
  tz_pos <- array(NA, c(n, p, mc))
  
  # MCMC iteration 
  pb <- progress_bar$new(total=mc)   # progress bar 
  for (item in 1:mc){
    # Q
    chol_solve_Si <- chol(solve_Si)
    solve_Si <- t(chol_solve_Si) %*% chol_solve_Si
    psi <- sqrt(diag(solve_Si))
    Q <- t(t(solve_Si) / psi)
    Q <- Q / psi
    
    # (s, r): local parameter
    for (i in 1:n){
      out <- (z[i,]==1)
      nonout <- !out
      l <- sum(z[i,])
      if( (p-l)>=1 ){
        s[i, nonout] <- sample(x=c(1, -1), size=p-l, replace=T, prob=rep(1/2, 2))
        vec <- rep(NA, p-l)
        for(h in 1:(p-l)) {
          vec[h] <- rejection(ga=ga, de=de_dash)
        }
        r[i, nonout] <- vec
        if( l>=1 ){
          u <- runif(l, min=0, max=1/(1+log(r[i,out]))^(1+ga))
          e_i <- Y[i,] 
          Precision_i <- e_i * t( e_i*t(solve_Si) )
          q_i <- as.list(rep(NA, 2*l))
          O_i <- matrix(0, l, l)
          for(h in 1:l){
            E_ih <- O_i
            E_ih[h, h] <- 1
            q_i[[h]] <- list(E_ih, rep(0, l), (-1)*exp( (1-1/(u[h])^(1/(1+ga)))*2 ))
            q_i[[l+h]] <- list((-1)*E_ih, rep(0, l), 1)
          }
          ze_tilde <- rnorm(l, mean=de_tilde*(1/t.[i, out]), sd=sqrt(de_tilde))
          setTimeLimit(cpu=5, elapsed=5)   # set limit time to avoid stacking
          value_tmg <- rtmg(n=1, M=(de_tilde+de_dash)*diag(l)+Precision_i[out, out, drop=F], 
                            r=ze_tilde+(-1)*colSums(Precision_i[nonout, out, drop=F]), 
                            initial=1/t.[i, out], f=NULL, g=NULL, q=q_i, burn.in=0)
          setTimeLimit(cpu=Inf, elapsed=Inf)
          value_tmg <- c(value_tmg)
          s[i, out] <- sign(value_tmg)
          r[i, out] <- 1/abs(value_tmg)
        }
      }else{
        u <- runif(l, min=0, max=1/(1+log(r[i, out]))^(1+ga))
        e_i <- Y[i,] 
        Precision_i <- e_i * t( e_i*t(solve_Si) )
        q_i <- as.list(rep(NA, 2*l))
        O_i <- matrix(0, l, l)
        for(h in 1:l){
          E_ih <- O_i
          E_ih[h, h] <- 1
          q_i[[h]] <- list(E_ih, rep(0, l), (-1)*exp( (1-1/(u[h])^(1/(1+ga)))*2 ))
          q_i[[l+h]] <- list((-1)*E_ih, rep(0, l), 1)
        }
        ze_tilde <- rnorm(l, mean=de_tilde*(1/t.[i, out]), sd=sqrt(de_tilde))
        setTimeLimit(cpu=5, elapsed=5)   # set limit time to avoid stacking
        value_tmg <- rtmg(n=1, M=(de_tilde+de_dash)*diag(l)+Precision_i[out, out, drop=F], 
                          r=ze_tilde, initial=1/t.[i, out], f=NULL, g=NULL, q=q_i, burn.in=0)
        setTimeLimit(cpu=Inf, elapsed=Inf)
        value_tmg <- c(value_tmg)
        s[i, out] <- sign(value_tmg)
        r[i, out] <- 1/abs(value_tmg)
      }
    }
    r[abs(r)>upper_bound] <- upper_bound*(sign(r))[abs(r)>upper_bound]
    if(positive){ 
      t. <- r
    }else{
      t. <- r*s 
    }
    tz <- (t.)^z
    
    # (s, z): local parameter 
    eigen_Q <- eigen(Q, symmetric=T, only.values=F)
    max_la <- eigen_Q$values[1]
    for (i in 1:n) {
      ze <- rnorm(p, sd=sqrt(de +max_la - eigen_Q$values), 
                  mean=(de + max_la - eigen_Q$values)*c(t(eigen_Q$vectors)%*%(psi*Y[i, ]/tz[i,])))
      C_i <- (phi/abs(r[i, ])) / (1-phi)
      al_p1 <- Y[i,] / (r[i, ]*1)^1
      al_p0 <- Y[i,] / (r[i, ]*1)^0
      al_m1 <- Y[i,] / (r[i, ]*(-1))^1
      al_m0 <- Y[i,] / (r[i, ]*(-1))^0
      log_den_p1 <- 1*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_p1 - 0.5*(max_la+de)*psi^2*al_p1^2
      log_den_p0 <- 0*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_p0 - 0.5*(max_la+de)*psi^2*al_p0^2
      log_den_m1 <- 1*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_m1 - 0.5*(max_la+de)*psi^2*al_m1^2
      log_den_m0 <- 0*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_m0 - 0.5*(max_la+de)*psi^2*al_m0^2
      log_den <- cbind(log_den_p1, log_den_p0, log_den_m1, log_den_m0)
      for (k in 1:p) {
        prob <- log_den[k, ]
        prob <- prob - max(prob)
        prob <- exp(prob)
        prob <- prob / sum(prob)
        number <- sample(x=1:4, size=1, replace=T, prob=prob)
        s[i, k] <- switch(number, first.=1, second=1, third.=-1, fourth=-1)
        z[i, k] <- switch(number, first.=1, second=0, third.=1, fourth=0)
      }
    }
    if(positive){ 
      t. <- r
    }else{
      t. <- r*s 
    }
    tz <- t.^z
    
    # phi (marginal outlier probability)
    sum_z <- sum(z)
    phi <- rbeta(1, shape1=sum_z+a, shape2=n*p-sum_z+b)
    
    # Si 
    cross_prod <- function(i){
      ga_i <- Y[i,]/tz[i,]
      return(outer(ga_i, ga_i))
    }
    Mat_outer <- sapply(1:n, cross_prod)
    Mat <- matrix(rowSums(Mat_outer), nrow=p, ncol=p)
    Si <- solve( rWishart(1, df=n+nu, Sigma=solve(Mat+Ga))[,,1] )
    solve_Si <- solve(Si)
    psi <- sqrt( diag(solve_Si) )
    Q <- t( t(solve_Si)/psi )
    Q <- Q / psi
    
    # save posterior samples 
    Si_pos[,,item] <- Si
    z_pos[,,item] <- z
    tz_pos[,,item] <- tz
    
    pb$tick()
  }
  
  # posterior samples 
  Si_pos <- Si_pos[,,-(1:bn)]
  z_pos <- z_pos[,,-(1:bn)]
  Om_pos <- array(NA, c(p, p, mc-bn))  # precision matrix
  for(item in 1:(mc-bn)){
    Om_pos[,,item] <- solve( Si_pos[,,item] )
  }
  
  # return
  Result <- list(Sigma=Si_pos, Z=z_pos, Om=Om_pos, tz=tz_pos)
  return(Result)
}





###--------------------------------------------###
###    Robust Gaussian Graphical Modeling      ###
###      with positive latent variables        ###
###   (SCMix: Sandwich Covariance Mixture)     ###
###--------------------------------------------###
SCMix_GM_positive <- function(Y, mc=2000, bn=1000){
  # preliminary
  upper_bound <- 10^8
  p <- ncol(Y)
  n <- nrow(Y)
  
  # hyper-parameters of priors 
  de <- 10^(-8)
  de_dash <- 0
  de_tilde <- 1
  ga <- 1
  a <- 0.05 
  b <- 1  
  nu <- p
  Ga <- (nu+p+1)*diag(p)
  
  # function (rejection)
  rejection <- function(ga, de){
    TF <- FALSE
    while(!TF){
      r <- rbeta(1, shape1=1, shape2=ga)
      r <- exp( r/(1-r) )
      TF <- ( log(runif(1)) <= (-de)/(2*r^2) )
    }
    return(r)
  }
  
  # initial values 
  r <- matrix(2, n, p)    # latent variable
  s <- matrix(1, n, p)    # latent variable
  z <- matrix(1, n, p)    # latent variable
  t. <- r*s
  tz <- t.^z
  phi <- 0.05    # outlier probability
  Si <- diag(p)
  solve_Si <- solve(Si)
  psi <- sqrt( diag(solve_Si) )
  Q <- t( t(solve_Si)/psi )
  Q <- Q/psi
  
  # objects for posterior samples
  Si_pos <- array(NA, c(p, p, mc))
  z_pos <- array(NA, c(n, p, mc))
  tz_pos <- array(NA, c(n, p, mc))
  
  # MCMC iteration 
  pb <- progress_bar$new(total=mc)   # progress bar 
  for (item in 1:mc){
    # Q
    chol_solve_Si <- chol(solve_Si)
    solve_Si <- t(chol_solve_Si) %*% chol_solve_Si
    psi <- sqrt(diag(solve_Si))
    Q <- t(t(solve_Si) / psi)
    Q <- Q / psi
    
    # Update of \tilde{r}_{i,k} (s_{i,k} is fixed to +1)
    for (i in 1:n){
      out <- (z[i,]==1)
      nonout <- !out
      l <- sum(z[i,])
      if( (p-l)>=1 ){
        vec <- rep(NA, p-l)
        for(h in 1:(p-l)) {
          vec[h] <- rejection(ga=ga, de=de_dash)
        }
        r[i, nonout] <- vec
        if( l>=1 ){
          u <- runif(l, min=0, max=1/(1+log(r[i,out]))^(1+ga))
          e_i <- Y[i,] 
          Precision_hat_i <- e_i * t( e_i*t(solve_Si) )
          q_i <- as.list(rep(NA, 2*l))
          O_i <- matrix(0, l, l)
          for(h in 1:l){
            E_ih <- O_i
            E_ih[h, h] <- 1
            q_i[[h]] <- list(E_ih, rep(0, l), (-1)*exp( (1-1/(u[h])^(1/(1+ga)))*2 ))
            q_i[[l+h]] <- list((-1)*E_ih, rep(0, l), 1)
          }
          ze_tilde <- rnorm(l, mean=de_tilde*(1/r[i, out]), sd=sqrt(de_tilde))
          setTimeLimit(cpu=5, elapsed=5)   # set limit time to avoid stacking
          value_tmg <- rtmg(n=1, M=(de_tilde+de_dash)*diag(l)+Precision_hat_i[out, out, drop=F], 
                            r=ze_tilde+(-1)*colSums(Precision_hat_i[nonout, out, drop=F]), 
                            initial=1/r[i, out], 
                            f = diag(l), g = rep(10^(-88), l), 
                            q=q_i, burn.in=0)
          setTimeLimit(cpu=Inf, elapsed=Inf)
          value_tmg <- c(value_tmg)
          value_tmg <- pmax(0, value_tmg)
          r[i, out] <- 1 / value_tmg
        }
      }else{
        u <- runif(l, min=0, max=1/(1+log(r[i,out]))^(1+ga))
        e_i <- Y[i,] 
        Precision_hat_i <- e_i * t( e_i*t(solve_Si) )
        q_i <- as.list(rep(NA, 2*l))
        O_i <- matrix(0, l, l)
        for(h in 1:l){
          E_ih <- O_i
          E_ih[h, h] <- 1
          q_i[[h]] <- list(E_ih, rep(0, l), (-1)*exp( (1-1/(u[h])^(1/(1+ga)))*2 ))
          q_i[[l+h]] <- list((-1)*E_ih, rep(0, l), 1)
        }
        ze_tilde <- rnorm(l, mean=de_tilde*(1/r[i, out]), sd=sqrt(de_tilde))
        setTimeLimit(cpu=5, elapsed=5)   # set limit time to avoid stacking
        value_tmg <- rtmg(n=1, M=(de_tilde+de_dash)*diag(l)+Precision_hat_i[out, out, drop=F], 
                          r=ze_tilde+(-1)*colSums(Precision_hat_i[nonout, out, drop=F]), 
                          initial=1/r[i, out], 
                          f = diag(l), g = rep(10^(-88), l), 
                          q=q_i, burn.in=0)
        setTimeLimit(cpu=Inf, elapsed=Inf)
        value_tmg <- c(value_tmg)
        value_tmg <- pmax(0, value_tmg)
        r[i, out] <- 1 / value_tmg
      }
    }
    r[r>upper_bound] <- upper_bound
    t. <- r
    tz <- (t.)^z
    
    # update of z (s_{ik} is fixed to +1)
    eigen_Q <- eigen(Q, symmetric=T, only.values=F)
    max_la <- eigen_Q$values[1]
    for (i in 1:n) {
      ze <- rnorm(p, sd=sqrt(de + max_la - eigen_Q$values), 
                  mean=(de + max_la - eigen_Q$values)*c(t(eigen_Q$vectors)%*%(psi*Y[i, ]/tz[i,])))
      C_i <- (phi/r[i, ]) / (1-phi)
      al_p1 <- Y[i,] / (r[i, ]*1)^1
      al_p0 <- Y[i,] / (r[i, ]*1)^0
      log_den_p1 <- 1*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_p1 - 0.5*(max_la+de)*psi^2*al_p1^2
      log_den_p0 <- 0*log(C_i) + c(eigen_Q$vectors%*%ze)*psi*al_p0 - 0.5*(max_la+de)*psi^2*al_p0^2
      log_den <- cbind(log_den_p1, log_den_p0)
      for (k in 1:p) {
        prob <- log_den[k, ]
        prob <- prob - max(prob)
        prob <- exp(prob)
        prob <- prob / sum(prob)
        z[i, k] <- ifelse(test=(runif(1)<=prob[1]), yes=1, no=0)
      }
    }
    t. <- r
    tz <- t.^z
    
    # phi (marginal outlier probability)
    sum_z <- sum(z)
    phi <- rbeta(1, shape1=sum_z+a, shape2=n*p-sum_z+b)
    
    # Si 
    cross_prod <- function(i){
      ga_i <- Y[i,]/tz[i,]
      return(outer(ga_i, ga_i))
    }
    Mat_outer <- sapply(1:n, cross_prod)
    Mat <- matrix(rowSums(Mat_outer), nrow=p, ncol=p)
    Si <- solve( rWishart(1, df=n+nu, Sigma=solve(Mat+Ga))[,,1] )
    solve_Si <- solve(Si)
    psi <- sqrt( diag(solve_Si) )
    Q <- t( t(solve_Si)/psi )
    Q <- Q / psi
    
    # save posterior samples 
    Si_pos[,,item] <- Si
    z_pos[,,item] <- z
    tz_pos[,,item] <- tz
    
    pb$tick()
  }
  
  # posterior samples 
  Si_pos <- Si_pos[,,-(1:bn)]
  z_pos <- z_pos[,,-(1:bn)]
  Om_pos <- array(NA, c(p, p, mc-bn))  # precision matrix
  for(item in 1:(mc-bn)){
    Om_pos[,,item] <- solve( Si_pos[,,item] )
  }
  
  # return
  Result <- list(Sigma=Si_pos, Z=z_pos, Om=Om_pos, tz=tz_pos)
  return(Result)
}





