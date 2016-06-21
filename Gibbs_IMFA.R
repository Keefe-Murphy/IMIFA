#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  gibbs.IMFA       <- function(Q, data, iters, N, P, G, mu.zero,
                               sigma.mu, burnin, thinning, mu, trunc.G,
                               psi.alpha, psi.beta, verbose, alpha.d1,
                               alpha.dk, sw, cluster, phi.nu, b0, b1, prop,
                               beta.d1, beta.dk, adapt, epsilon, ...) {
        
  # Define & initialise variables
    n.iters        <- round(max(iters), -1)
    n.store        <- length(iters)
    Gs             <- seq_len(G)
    Ts             <- seq_len(trunc.G)
    Ps             <- seq_len(P)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Ts)
    iternames      <- paste0("Iteration", seq_len(n.store))
    Q0             <- Q > 0
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, trunc.G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store      <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, trunc.G, n.store))
      dimnames(load.store) <- list(varnames, if(Q0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store    <- array(0, dim=c(P, trunc.G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store     <- matrix(0, nr=trunc.G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store        <- matrix(0, nr=N, nc=n.store)
    ll.store       <- rep(0, n.store)
    G.store        <- rep(0, n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    
    mu.sigma       <- 1/sigma.mu
    l.sigma        <- 1/sigma.l 
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gs)
    pi.alpha       <- cluster$pi.alpha + N
    pi.prop        <- cbind(cluster$pi.prop, sim.stick(pi.alpha=cluster$pi.alpha, nn=rep(0, trunc.G))[,-Gs, drop=F])
    mu             <- cbind(mu, do.call(cbind, lapply(seq_len(trunc.G - G), function(g) sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero))))
    f              <- sim.f.p(N=N, Q=Q)
    lmat           <- lapply(Ts, function(t) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F))
    psi.inv        <- do.call(cbind, lapply(Ts, function(t) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    if(Q0) {
      for(g in Gs)   {
        fact       <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } else                warning(paste0("Parameters of group ", g, " initialised by simulation from priors, not factanal: G=", G, ", Q=", Q), call.=F)
      }
    } else {
      psi.inv[,Gs]         <- do.call(cbind, lapply(Gs, function(g) if(pi.prop[,g] > 0) 1/apply(data[z == g,, drop=F], 2, var) else psi.inv[,g]))
    }
    l.sigma        <- l.sigma * diag(Q)
    lmat           <- array(unlist(lmat, use.names=F), dim=c(P, Q, trunc.G))
    index          <- order(pi.prop, decreasing=TRUE)
    pi.prop        <- pi.prop[,index, drop=F]
    mu             <- mu[,index, drop=FALSE]
    lmat           <- lmat[,,index, drop=FALSE]
    psi.inv        <- psi.inv[,index, drop=FALSE]
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, G=G, pi.prop=pi.prop, Sigma=lapply(Gs,
                                  function(g) tcrossprod(as.matrix(lmat[,,g])) + diag(1/psi.inv[,g])))$log.likes)
      G.store[1]           <- G
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
    
    # Slice Sampler
      nn           <- tabulate(z, nbins=trunc.G)
      pi.i         <- pi.prop[,z]
      u.slice      <- runif(N, 0, pi.i)
      Au           <- unlist(lapply(seq_along(u.slice), function(u) sum(u.slice[u] < pi.prop)))
      G            <- max(Au)
      Gs           <- seq_len(G)
      nn0          <- nn[Gs] > 0
      z.ind        <- lapply(Gs, function(g) z == g)
      slice.ind    <- do.call(cbind, lapply(Gs, function(g) u.slice < pi.prop[,g]))
      
    # Means
      sum.data     <- lapply(Gs, function(g) colSums(data[z.ind[[g]],, drop=F]))
      sum.f        <- lapply(Gs, function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu[,Gs]      <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], 
                              sum.f=sum.f[[g]], lmat=as.matrix(lmat[,,g]), mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero)))
      
    # Scores & Loadings
      c.data       <- lapply(Gs, function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(Q0)   {
        f          <- do.call(rbind, lapply(Gs, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=as.matrix(lmat[,,g]), 
                              c.data=c.data[[g]], psi.inv=psi.inv[,g], Q=Q)))[obsnames,, drop=F]
        FtF        <- lapply(Gs, function(g) if(nn0[g]) crossprod(f[z.ind[[g]],, drop=F]))
        lmat[,,Gs] <- array(unlist(lapply(Gs, function(g) if(nn0[g]) matrix(unlist(lapply(Ps, function(j) sim.load(l.sigma=l.sigma, Q=Q, P=P, c.data=c.data[[g]][,j],  f=f[z.ind[[g]],, drop=F], 
                            psi.inv=psi.inv[,g][j], FtF=FtF[[g]], shrink=F)), use.names=F), nr=P, byrow=T) else sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F)), use.names=F), dim=c(P, Q, G))
      }
                    
    # Uniquenesses
      psi.inv[,Gs] <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.psi.i(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, 
                              P=P, f=f[z.ind[[g]],,drop=F], lmat=as.matrix(lmat[,,g])) else sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
      
    # Mixing Proportions
      pi.prop      <- sim.stick(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi          <- 1/psi.inv
      Sigma        <- lapply(Gs, function(g) tcrossprod(as.matrix(lmat[,,g])) + diag(psi[,g]))
      z.res        <- sim.z(data=data, mu=mu, Sigma=Sigma, G=G, pi.prop=pi.prop, slice.ind=slice.ind)
      z            <- z.res$z
    
    # Label Switching
      switch.lab   <- lab.switch(z.new=z, z.old=z.temp, Gs=Gs)
      z            <- switch.lab$z
      z.perm       <- switch.lab$z.perm
      perm         <- identical(as.integer(z.perm), Gs)
      if(!perm) {
        mu         <- mu[,z.perm, drop=F]
        lmat       <- lmat[,,z.perm, drop=F]
        psi.inv    <- psi.inv[,z.perm, drop=F]
        pi.prop    <- pi.prop[,z.perm, drop=F]
      }
    
      if(is.element(iter, iters))   {
        new.it     <- which(iters == iter)
        log.like   <- sum(z.res$log.likes)
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu 
        if(all(sw["f.sw"], Q0))    f.store[,,new.it]       <- f
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
                                   z.store[,new.it]        <- z 
                                   ll.store[new.it]        <- log.like
                                   G.store[new.it]         <- sum(nn0)
      }
    }
    returns        <- list(mu       = if(sw["mu.sw"])         mu.store,
                           f        = if(all(sw["f.sw"], Q0)) f.store, 
                           load     = if(all(sw["l.sw"], Q0)) load.store, 
                           psi      = if(sw["psi.sw"])        psi.store,
                           pi.prop  = if(sw["pi.sw"])         pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store)
    return(returns)
  }