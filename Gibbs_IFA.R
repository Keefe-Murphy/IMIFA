###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Gibbs Sampler Function
  gibbs.IFA      <- function(Q = NULL, data = NULL, iters = NULL, 
                             N = NULL, P = NULL, sigma.mu = NULL, 
                             psi.alpha = NULL, psi.beta = NULL, 
                             burnin = NULL, thinning = NULL, verbose = NULL, 
                             sw = NULL, phi.nu = NULL, alpha.d1 = NULL, 
                             alpha.d2 = NULL, adapt = NULL, b0 = NULL, 
                             b1 = NULL, prop = NULL, epsilon = NULL, ...) {    
    
  # Define & initialise variables
    n.iters      <- round(max(iters), -1)
    n.store      <- length(iters)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    iternames    <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store   <- matrix(0, nr=P, nc=n.store)
      dimnames(mu.store)   <- list(varnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store  <- matrix(0, nr=P, nc=n.store)
      dimnames(psi.store)  <- list(varnames, iternames)
    }
    post.mu      <- setNames(rep(0, P), varnames)
    post.psi     <- setNames(rep(0, P), varnames)
    cov.emp      <- cov(data)
    cov.est      <- matrix(0, nr=P, nc=P)
    Q.star       <- Q
    Q.store      <- setNames(rep(0, n.store), iternames)
    dimnames(cov.emp)      <- list(varnames, varnames)
    dimnames(cov.est)      <- dimnames(cov.emp)
    
    mu.sigma     <- 1/sigma.mu
    mu           <- sim.mu.p(P=P, sigma.mu=sigma.mu)  
    f            <- sim.f.p(Q=Q, N=N)
    psi.inv      <- sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    phi          <- sim.phi.p(Q=Q, P=P, phi.nu=phi.nu)
    delta        <- sim.delta.p(Q=Q, alpha.d1=alpha.d1, alpha.d2=alpha.d2)
    tau          <- cumprod(delta)
    lmat         <- sim.load.p(Q=Q, phi=phi, tau=tau, P=P)
    sum.data     <- colSums(data)
    if(burnin     < 1) {
      mu.store[,1]         <- mu
      f.store[,,1]         <- f
      load.store[,,1]      <- lmat
      psi.store[,1]        <- 1/psi.inv
    }
  
  # Iterate
    for(iter in seq_len(max(iters))[-1]) { 
      if(verbose) {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
    # Means
      sum.f      <- colSums(f)
      mu         <- sim.mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv,
                           sum.data=sum.data, sum.f=sum.f, lmat=lmat)
      
    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q  > 0) {
        f        <- sim.score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data)
        FtF      <- crossprod(f)
        lmat     <- sim.load(Q=Q, c.data=c.data, P=P, f=f, phi=phi,
                             tau=tau, psi.inv=psi.inv, FtF=FtF)
      } else {
        f        <- matrix(, nr=N, nc=0)
        lmat     <- matrix(, nr=P, nc=0)
      }          
    
    # Uniquenesses
      psi.inv    <- sim.psi.i(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                              c.data=c.data, f=f, lmat=lmat)
    
    # Local Shrinkage
      load.2     <- lmat * lmat
      phi        <- sim.phi(Q=Q, P=P, phi.nu=phi.nu, tau=tau, load.2=load.2)
          
    # Global Shrinkage
      sum.term   <- diag(crossprod(phi, load.2))
    if(Q  > 0) {
      delta[1]   <- sim.delta1(Q=Q, P=P, alpha.d1=alpha.d1, delta=delta,
                               tau=tau, sum.term=sum.term)
      tau        <- cumprod(delta)  
    } 
    if(Q >= 2) {
      for(k in seq_len(Q)[-1]) { 
        delta[k] <- sim.deltak(Q=Q, P=P, k=k, alpha.d2=alpha.d2,
                               delta=delta, tau=tau, sum.term=sum.term)
        tau      <- cumprod(delta)      
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
     if(is.element(iter, iters)) {
        new.it   <- which(iters == iter)  
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        Sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + Sigma/n.store
        if(sw["mu.sw"])             mu.store[,new.it]              <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,seq_len(Q),new.it]    <- f
        if(all(sw["l.sw"], Q > 0))  load.store[,seq_len(Q),new.it] <- lmat
        if(sw["psi.sw"])            psi.store[,new.it]             <- psi
                                    Q.store[new.it]                <- Q
      }
    }
    returns      <- list(mu       = if(sw["mu.sw"])  mu.store,
                         f        = if(sw["f.sw"])   as.simple_sparse_array(f.store), 
                         load     = if(sw["l.sw"])   as.simple_sparse_array(load.store), 
                         psi      = if(sw["psi.sw"]) psi.store,
                         post.mu  = post.mu,
                         post.psi = post.psi,
                         cov.emp  = cov.emp,
                         cov.est  = cov.est,
                         Q.store  = Q.store)
    return(returns)
  }