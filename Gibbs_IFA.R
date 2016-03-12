###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Gibbs Sampler Function
  gibbs.IFA      <- function(Q, data, n.iters, N, P, 
                             sigma.mu, psi.alpha, psi.beta, 
                             burnin, thinning, n.store, verbose, 
                             sw, phi.nu, alpha.d1, alpha.d2, 
                             adapt, b0, b1, prop, epsilon, ...) {    
    
  # Define & initialise variables
    cnames       <- colnames(data)
    rnames       <- rownames(data)
    facnames     <- paste0("Factor ", 1:Q)
    iternames    <- paste0("Iteration", 1:n.store)
    if(sw["mu.sw"]) {
      mu.store   <- matrix(0, nr=P, nc=n.store)
      dimnames(mu.store)     <- list(cnames, iternames)
    }
    if(sw["f.sw"])  {
      f.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)      <- list(rnames, facnames, iternames)
    }
    if(sw["l.sw"])  {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store)   <- list(cnames, facnames, iternames)
    }
    if(sw["p.sw"])  {
      psi.store  <- matrix(0, nr=P, nc=n.store)
      dimnames(psi.store)    <- list(cnames, iternames)
    }
    post.mu      <- setNames(rep(0, P), cnames)
    post.psi     <- setNames(rep(0, P), cnames)
    post.Sigma   <- matrix(0, nr=P, nc=P)
    cov.emp      <- cov(data)
    Q.star       <- Q
    Q.store      <- setNames(rep(0, n.store), iternames)
    dimnames(post.Sigma)     <- list(cnames, cnames)
    dimnames(cov.emp)        <- dimnames(post.Sigma)
    
    sigma.mu     <- 1/sigma.mu
    mu           <- sim.mu.p(sigma.mu, P)  
    f            <- sim.f.p(Q, N)
    psi.inv      <- sim.pi.p(P, psi.alpha, psi.beta)
    phi          <- sim.p.p(Q, P, phi.nu)
    delta        <- sim.d.p(Q, alpha.d1, alpha.d2)
    tau          <- cumprod(delta)
    lmat         <- matrix(0, nr=P, nc=Q)
    for(j in 1:P) {
      D.load     <- phi[j,] * tau
      lmat[j,]   <- sim.l.p(D.load, Q)
    }
    sum.data     <- colSums(data)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(verbose) {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
    # Means
      sum.f      <- colSums(f)
      mu         <- sim.mu(N, P, sigma.mu, psi.inv, sum.data, sum.f, lmat)
      
    # Scores
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q  > 0) {
        f        <- sim.scores(N, Q, lmat, psi.inv, c.data)
      } else {
        f        <- matrix(, nr=N, nc=0)
      }          
    
    # Loadings
      FtF        <- crossprod(f)
      if(Q  > 0) {
        for (j in 1:P) {
          psi.inv.j <- psi.inv[j]
          c.data.j  <- c.data[,j]
          D.load    <- phi[j,] * tau * diag(Q)
          lmat[j,]  <- sim.load(D.load, Q, c.data.j, f, psi.inv.j, FtF)
        } 
      } else {
        lmat     <- matrix(, nr=P, nc=0)
      }
    
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N, P, psi.alpha, psi.beta, c.data, f, lmat)
    
    # Local Shrinkage
      load.2     <- lmat * lmat
      phi        <- sim.phi(Q, P, phi.nu, tau, load.2)
          
    # Global Shrinkage
      sum.term   <- diag(crossprod(phi, load.2))
      if(Q  > 0) {
        delta[1]    <- sim.delta1(Q, P, alpha.d1, delta, tau, sum.term)
        tau         <- cumprod(delta)  
      } 
      if(Q >= 2) {
        for(k in 2:Q) { 
          delta[k]  <- sim.deltak(Q, P, k, alpha.d2, delta, tau, sum.term)
          tau       <- cumprod(delta)      
        }
      }
    
    # Adaptation  
      if(all(adapt, iter > burnin)) {      
        prob     <- 1/exp(b0 + b1 * pmax(iter - burnin, 0))
        unif     <- runif(n=1, min=0, max=1)
        if(Q > 0) {
          lind   <- colSums(abs(lmat) < epsilon) / P
        } else {
          lind   <- 0
        }
        colvec   <- lind >= prop
        numred   <- sum(colvec)
        
        if(unif   < prob) { # check whether to adapt or not
          if(numred == 0) { # simulate extra columns from priors
            Q       <- Q + 1
            f       <- cbind(f, rnorm(n=N, mean=0, sd=1))         
            phi     <- cbind(phi, rgamma(n=P, shape=phi.nu/2, rate=phi.nu/2))
            delta   <- c(delta, rgamma(n=1, shape=alpha.d2, rate=1))
            tau     <- cumprod(delta)
            lmat    <- cbind(lmat, rnorm(n=P, mean=0, sd=sqrt(1/phi[,Q] * 1/tau[Q])))
          } else          { # remove redundant columns
            nonred  <- which(colvec == 0)
            Q       <- Q - numred
            f       <- f[,nonred, drop=F]
            phi     <- phi[,nonred, drop=F]
            delta   <- delta[nonred]
            tau     <- cumprod(delta)
            lmat    <- lmat[,nonred, drop=F]
          }
        }
      } 
    
    if(Q > Q.star)  stop(paste0("Q cannot exceed initial number of loadings columns: try increasing Q.star from ", Q.star))
      if(all(iter > burnin, iter %% thinning == 0)) {
        new.iter <- ceiling((iter - burnin)/thinning)
        psi      <- 1/psi.inv
        if(sw["mu.sw"])             mu.store[,new.iter]       <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,1:Q,new.iter]    <- f
        if(all(sw["l.sw"], Q > 0))  load.store[,1:Q,new.iter] <- lmat
        if(sw["p.sw"])              psi.store[,new.iter]      <- psi
        post.mu     <- post.mu + mu/n.store
        post.psi    <- post.psi + psi/n.store
        Sigma       <- tcrossprod(lmat) + diag(psi)
        post.Sigma  <- post.Sigma + Sigma/n.store
        Q.store[new.iter]    <- Q
      }
    }
    returns      <- list(mu   = if(sw["mu.sw"]) mu.store,
                         f    = if(sw["f.sw"])  as.simple_sparse_array(f.store), 
                         load = if(sw["l.sw"])  as.simple_sparse_array(load.store), 
                         psi  = if(sw["p.sw"])  psi.store,
                         cov.mat    = cov.emp,
                         post.mu    = post.mu,
                         post.psi   = post.psi,
                         post.Sigma = post.Sigma,
                         Q.store    = Q.store)
    return(returns[!sapply(returns, is.null)])
  }