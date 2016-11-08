################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Gibbs Sampler Function
  gibbs.FA       <- function(Q, data, iters, N, P, sigma.mu, mu,
                             mu.zero, psi.alpha, psi.beta, burnin,
                             thinning, verbose, sw, sigma.l, ...) {
        
  # Define & initialise variables
    start.time   <- proc.time()
    total        <- max(iters)
    if(verbose)     pb     <- txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    iternames    <- paste0("Iteration", seq_len(n.store))
    Q0           <- Q  > 0
    Q1           <- Q == 1
    dimnames(data)         <- NULL
    if(sw["mu.sw"])  {
      mu.store   <- matrix(0, nr=P, nc=n.store)
      dimnames(mu.store)   <- list(varnames, iternames)
    }
    if(sw["s.sw"])   {
      eta.store  <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)  <- list(obsnames, if(Q0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store) <- list(varnames, if(Q0) facnames, iternames)
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
    dimnames(cov.emp)      <- list(varnames, varnames)
    dimnames(cov.est)      <- dimnames(cov.emp)
    
    mu.sigma     <- 1/sigma.mu
    eta          <- sim.eta.p(Q=Q, N=N)
    lmat         <- sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE)
    psi.inv      <- sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    if(all(Q0, Q  < P - sqrt(P + Q), N > P)) {
      fact       <- try(factanal(data, factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
      if(!inherits(fact, "try-error")) {
        eta      <- fact$scores
        lmat     <- fact$loadings
        psi.inv  <- 1/fact$uniquenesses
      }
    } else {
      psi.tmp    <- psi.inv
      psi.inv    <- 1/apply(data, 2, var)
      inf.ind    <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma      <- diag(1/sigma.l, Q)
    sum.data     <- mu * N
    if(burnin     < 1)    {
      mu.store[,1]         <- mu
      eta.store[,,1]       <- eta
      load.store[,,1]      <- lmat
      psi.store[,1]        <- 1/psi.inv
      ll.store[1]          <- sum(dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + diag(1/psi.inv), log=TRUE))
    }
    init.time    <- proc.time() - start.time
  
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose && iter    < burnin) setTxtProgressBar(pb, iter)
    
    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q0) {
        eta      <- sim.score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(unlist(lapply(Pseq, function(j) sim.load(l.sigma=l.sigma, Q=Q, eta=eta, c.data=c.data[,j], P=P, 
                           Q1=Q1, psi.inv=psi.inv[j], EtE=crossprod(eta), shrink=FALSE)), use.names=FALSE), nr=P, byrow=TRUE)
      }
      
    # Means
      mu[]       <- sim.mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums(eta), lmat=lmat, mu.zero=mu.zero)
                      
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, c.data=c.data, eta=eta, lmat=lmat)
    
      if(is.element(iter, iters)) {
        if(verbose) setTxtProgressBar(pb, iter)
        new.it   <- which(iters == iter)  
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + sigma/n.store
        if(sw["mu.sw"])          mu.store[,new.it]    <- mu  
        if(all(sw["s.sw"], Q0))  eta.store[,,new.it]  <- eta
        if(all(sw["l.sw"], Q0))  load.store[,,new.it] <- lmat
        if(sw["psi.sw"])         psi.store[,new.it]   <- psi
                                 ll.store[new.it]     <- sum(dmvn(X=data, mu=mu, sigma=sigma, log=TRUE))
      }  
    }
    close(pb)
    returns   <- list(mu       = if(sw["mu.sw"])         mu.store,
                      eta      = if(all(sw["s.sw"], Q0)) eta.store, 
                      load     = if(all(sw["l.sw"], Q0)) load.store, 
                      psi      = if(sw["psi.sw"])        psi.store,
                      post.mu  = post.mu,
                      post.psi = post.psi,
                      cov.emp  = cov.emp,
                      cov.est  = cov.est,
                      ll.store = ll.store,
                      time     = init.time)
    attr(returns, "K")        <- P * Q - 0.5 * Q * (Q - 1) + 2 * P
    return(returns)
  }