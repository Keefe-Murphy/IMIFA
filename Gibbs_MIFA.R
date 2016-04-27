################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MIFA      <- function(Q = NULL, data = NULL, iters = NULL, 
                              N = NULL, P = NULL, G = NULL, sigma.mu = NULL, 
                              burnin = NULL, thinning = NULL, psi.alpha = NULL, 
                              psi.beta = NULL, sw = NULL, verbose = NULL, clust = NULL, 
                              mu0g = NULL, phi.nu = NULL, alpha.d1 = NULL, alpha.d2 = NULL, 
                              adapt = NULL, b0 = NULL, b1 = NULL, prop = NULL, epsilon = NULL, ...) {
        
  # Define & initialise variables
    n.iters       <- round(max(iters), -1)
    n.store       <- length(iters)
    obsnames      <- rownames(data)
    varnames      <- colnames(data)
    Gseq          <- seq_len(G)
    facnames      <- paste0("Factor ", seq_len(Q))
    gnames        <- paste0("Group ", Gseq)
    iternames     <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store    <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store     <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store  <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store   <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store    <- matrix(0, nr=G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store       <- matrix(0, nr=N, nc=n.store)
    ll.store      <- rep(0, n.store)
    fin.ll        <- T
    Q.star        <- Q
    Q.store       <- matrix(0, nr=G, nc=n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    dimnames(Q.store)      <- list(gnames, iternames)
    
    mu.sigma      <- 1/sigma.mu
    l.sigma       <- 1/sigma.l 
    pi.alpha      <- clust$pi.alpha
    z             <- clust$z
    pi.prop       <- t(prop.table(tabulate(z, nbins=G)))
    mu            <- do.call(cbind, lapply(Gseq, function(g) colMeans(data[z == g,, drop=F])))
    f             <- sim.f.p(N=N, Q=Q)
    fg            <- list()
    phi           <- lapply(Gseq, function(g) sim.phi.p(Q=Q, P=P, phi.nu=phi.nu))
    delta         <- lapply(Gseq, function(g) sim.delta.p(Q=Q, alpha.d1=alpha.d1, alpha.d2=alpha.d2))
    tau           <- lapply(delta, cumprod)
    lmat          <- lapply(Gseq, function(g) sim.load.p(Q=Q, phi=phi[[g]], tau=tau[[g]], P=P))
    psi.inv       <- do.call(cbind, lapply(Gseq, function(g) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    if(mu0g)  {
      mu.zero     <- mu
    } else    {
      mu.zero     <- do.call(cbind, lapply(Gseq, function(g) colMeans(data)))
    }
    if(Q > 0) {
      for(g in Gseq) {
        fact      <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } else                warning(paste0("Parameters of group ", g, " initialised by simulation from priors, not factanal: G=", G, ", Q=", Q), call.=F)
      }
    } else {
      psi.inv     <- 1/do.call(cbind, lapply(Gseq, function(g) apply(data[z == g,, drop=F], 2, var)))
    }
    l.sigma       <- l.sigma * diag(Q)
    Qs            <- rep(Q, G)
    if(burnin      < 1)     {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      for(g in Gseq) {
        load.store[,,g,1]  <- lmat[[g]]
      }
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
    }
    
  # Iterate
    for(iter in seq_len(max(iters))[-1]) { 
      if(verbose)  {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      nn          <- tabulate(z, nbins=G)
      z.ind       <- lapply(Gseq, function(g) z == g)
      
    # Means
      sum.data    <- lapply(Gseq, function(g) colSums(data[z.ind[[g]],,drop=F]))
      sum.f       <- lapply(Gseq, function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu          <- do.call(cbind, lapply(Gseq, function(g) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], 
                             P=P, sum.data=sum.data[[g]], sum.f=sum.f[[g]], lmat=lmat[[g]], mu.zero=mu.zero[,g])))
    
    # Scores & Loadings
      c.data      <- lapply(Gseq, function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(all(Qs   == 0))  {
        f         <- matrix(, nr=N, nc=0)
        lmat      <- lapply(Gseq, function(g) matrix(, nr=P, nc=0))
      } else {
        for(g in Gseq)    {
          Qg      <- Qs[g]
          c.datg  <- c.data[[g]]
          psi.ig  <- psi.inv[,g]
          if(Qg    > 0)   {
            fgg            <- sim.score(N=nn[g], lmat=lmat[[g]], Q=Qg, 
                                        c.data=c.datg, psi.inv=psi.ig)
            FtF            <- crossprod(fgg)
            lmat[[g]]      <- sim.load(l.sigma=l.sigma, Q=Qg, c.data=c.datg, P=P, f=fgg,
                                       psi.inv=psi.ig, FtF=FtF, phi=phi[[g]], tau=tau[[g]])
            fg[[g]]        <- fgg
          } else {
            fg[[g]]        <- matrix(, nr=nn[g], nc=0)
            lmat[[g]]      <- matrix(, nr=P, nc=0)
          }
        }
        f         <- do.call(rbind, fg)[obsnames,, drop=F]
      }
                  
    # Uniquenesses
      psi.inv     <- do.call(cbind, lapply(Gseq, function(g) sim.psi.i(N=nn[g], P=P, psi.alpha=psi.alpha, 
                             psi.beta=psi.beta, c.data=c.data[[g]], f=f[z.ind[[g]],,drop=F], lmat=lmat[[g]])))
    
    # Local Shrinkage
      load.2      <- lapply(lmat, function(x) x * x)
      phi         <- lapply(Gseq, function(g) sim.phi(Q=Qs[g], P=P, phi.nu=phi.nu, tau=tau[[g]], load.2=load.2[[g]]))
    
    # Global Shrinkage
      sum.term    <- lapply(Gseq, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gseq)    {
        Qg        <- Qs[g]
        sum.termg <- sum.term[[g]]
        if(Qg      > 0) {
          delta[[g]][1]    <- sim.delta1(Q=Qg, alpha.d1=alpha.d1, delta=delta[[g]], 
                                         P=P, tau=tau[[g]], sum.term=sum.termg)
          tau[[g]]         <- cumprod(delta[[g]])
        }
        if(Qg     >= 2) {
          for(k in seq_len(Qg)[-1]) { 
            delta[[g]][k]  <- sim.deltak(Q=Qg, alpha.d2=alpha.d2, delta=delta[[g]], 
                                         P=P, k=k, tau=tau[[g]], sum.term=sum.termg)
            tau[[g]]       <- cumprod(delta[[g]])
          }
        }
      }
    
    # Mixing Proportions
      pi.prop    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi        <- 1/psi.inv
      Sigma      <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      z.res      <- sim.z(data=data, mu=mu, Sigma=Sigma, G=G, pi.prop=pi.prop)
      z          <- z.res$z
      
    if(any(Qs > Q.star))      stop(paste0("Q cannot exceed initial number of loadings columns: try increasing Q.star from ", Q.star))
      if(is.element(iter, iters))  {
        new.it   <- which(iters == iter)
        log.like <- sum(z.res$log.likes)
        if(all(!is.finite(log.like),
           isTRUE(fin.ll))) { warning("Infinite likelihood: model-selection criteria may not be obtainable", call.=F)
          fin.ll <- F
        }
        if(sw["mu.sw"])             mu.store[,,new.it]     <- mu  
        if(all(sw["f.sw"], 
           any(Qs > 0)))            f.store[,,new.it]      <- f
        if(sw["l.sw"]) {
          for(g in Gseq)    {
            if(Qs[g] > 0)   {       load.store[,,g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])            psi.store[,,new.it]    <- psi
        if(sw["pi.sw"])             pi.store[,new.it]      <- pi.prop
                                    z.store[,new.it]       <- z 
                                    ll.store[new.it]       <- log.like
      }  
    }
    returns   <- list(mu       = if(sw["mu.sw"])    mu.store,
                      f        = if(all(sw["f.sw"]) as.simple_sparse_array(f.store), 
                      load     = if(all(sw["l.sw"]) as.simple_sparse_array(load.store), 
                      psi      = if(sw["psi.sw"])   psi.store,
                      pi.prop  = if(sw["pi.sw"])    pi.store,
                      z        = z.store,
                      ll.store = ll.store)
    return(returns)
  }