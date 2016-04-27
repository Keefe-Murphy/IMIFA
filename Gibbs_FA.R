################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Gibbs Sampler Function
  gibbs.FA       <- function(Q = NULL, data = NULL, iters = NULL,
                             N = NULL, P = NULL, sigma.mu = NULL,
                             psi.alpha = NULL, psi.beta = NULL,
                             burnin = NULL, thinning = NULL, 
                             verbose = NULL, sw = NULL, sigma.l = NULL, ...) {
        
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
    ll.store     <- rep(0, n.store)
    fin.ll       <- T
    dimnames(cov.emp)      <- list(varnames, varnames)
    dimnames(cov.est)      <- dimnames(cov.emp)
    
    mu.sigma     <- 1/sigma.mu
    l.sigma      <- 1/sigma.l
    mu           <- sim.mu.p(P=P, sigma.mu=sigma.mu)  
    f            <- sim.f.p(Q=Q, N=N)
    lmat         <- sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F)
    psi.inv      <- sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    l.sigma      <- l.sigma * diag(Q)
    sum.data     <- colSums(data)
    if(burnin     < 1)    {
      mu.store[,1]         <- mu
      f.store[,,1]         <- f
      load.store[,,1]      <- lmat
      psi.store[,1]        <- 1/psi.inv
      ll.store[1]          <- sum(mvdnorm(data=data, mu=mu, Sigma=tcrossprod(lmat) + diag(1/psi.inv), log.d=T))
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
      if(Q > 0) {
        f        <- sim.score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data)
        FtF      <- crossprod(f)
        lmat     <- sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data, f=f, 
                             P=P, psi.inv=psi.inv, FtF=FtF, shrink=F)
      }
                      
    # Uniquenesses
      psi.inv    <- sim.psi.i(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                              c.data=c.data, f=f, lmat=lmat)
    
      if(is.element(iter, iters)) {
        new.it   <- which(iters == iter)  
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        Sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + Sigma/n.store
        log.like <- sum(mvdnorm(data=data, mu=mu, Sigma=Sigma, log.d=T))
        if(all(!is.finite(log.like),
           isTRUE(fin.ll))) { warning("Infinite likelihood: model-selection criteria may not be obtainable", call.=F)
          fin.ll <- F
        }
        if(sw["mu.sw"])             mu.store[,new.it]    <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,new.it]    <- f
        if(all(sw["l.sw"], Q > 0))  load.store[,,new.it] <- lmat
        if(sw["psi.sw"])            psi.store[,new.it]   <- psi
                                    ll.store[new.it]     <- log.like
      }  
    }
    returns   <- list(mu       = if(sw["mu.sw"])            mu.store,
                      f        = if(all(sw["f.sw"], Q > 0)) f.store, 
                      load     = if(all(sw["l.sw"], Q > 0)) load.store, 
                      psi      = if(sw["psi.sw"])           psi.store,
                      post.mu  = post.mu,
                      post.psi = post.psi,
                      cov.emp  = cov.emp,
                      cov.est  = cov.est,
                      ll.store = ll.store)
    attr(returns, "K")        <- P * Q - 0.5 * Q * (Q - 1) + 2 * P
    return(returns)
  }