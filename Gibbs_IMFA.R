#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  gibbs.IMFA       <- function(Q, data, iters, N, P, G, mu.zero, pp, sigma.l, MH.step,
                               sigma.mu, burnin, thinning, mu, trunc.G, gen.slice, MH.lower,
                               psi.alpha, psi.beta, verbose, sw, cluster, MH.upper, ...) {
        
  # Define & initialise variables
    n.iters        <- round(max(iters), -1)
    n.store        <- length(iters)
    Gs             <- seq_len(G)
    Ts             <- seq_len(trunc.G)
    Ps             <- seq_len(P)
    Ns             <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Ts)
    iternames      <- paste0("Iteration", seq_len(n.store))
    Q0             <- Q > 0
    Q0s            <- rep(Q0, trunc.G)
    Q1             <- Q > 1
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
    if(MH.step)   {
      rate         <- rep(0, n.store)
      alpha.store  <- rep(0, n.store)
    }
    mu.sigma       <- 1/sigma.mu
    l.sigma        <- 1/sigma.l 
    z              <- cluster$z
    pi.alpha       <- cluster$pi.alpha
    pi.prop        <- cbind(cluster$pi.prop, sim.pi(pi.alpha=pi.alpha, nn=rep(0, trunc.G, inf.G=T))[,-Gs, drop=F])
    nn             <- tabulate(z, nbins=trunc.G)
    mu             <- cbind(mu, do.call(cbind, lapply(seq_len(trunc.G - G), function(g) sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero))))
    f              <- sim.f.p(N=N, Q=Q)
    lmat           <- lapply(Ts, function(t) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F))
    psi.inv        <- do.call(cbind, lapply(Ts, function(t) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    if(Q0) {
      for(g in which(nn > 2.5 * Q))      {
        fact       <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } 
      }
    } else {
      psi.tmp      <- psi.inv
      psi.inv      <- do.call(cbind, lapply(Ts, function(t) if(nn[g] > 1) 1/apply(data[z == g,, drop=F], 2, var) else psi.tmp[,g]))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[is.infinite(psi.inv)]
    }
    l.sigma        <- l.sigma * diag(Q)
    lmat           <- array(unlist(lmat, use.names=F), dim=c(P, Q, trunc.G))
    index          <- order(pi.prop, decreasing=T)
    pi.prop        <- pi.prop[,index, drop=F]
    mu             <- mu[,index, drop=F]
    lmat           <- lmat[,,index, drop=F]
    psi.inv        <- psi.inv[,index, drop=F]
    csi            <- pp * (1 - pp)^(Ts - 1)
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gs, N=N, pi.prop=pi.prop, Sigma=lapply(Gs,
                                  function(g) tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])), Q0=Q0s)$log.likes)
      G.store[1]           <- sum(nn > 0)
      if(MH.step)  {
        rate[1]            <- 0
        alpha.store[1]     <- pi.alpha 
      }
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
    
    # Mixing Proportions
      weights      <- sim.pi(pi.alpha=pi.alpha, nn=nn, inf.G=T, len=trunc.G)
      pi.prop      <- weights$pi.prop
      if(MH.step)     Vs   <- weights$Vs
      
    # Slice Sampler
      if(!gen.slice) {
        index      <- order(pi.prop, decreasing=T)
        pi.prop    <- pi.prop[,index, drop=F]
        mu         <- mu[,index, drop=F]
        lmat       <- lmat[,,index, drop=F]
        psi.inv    <- psi.inv[,index, drop=F]
        csi        <- pi.prop
      }
      u.slice      <- runif(N, 0, csi[z])
      G            <- max(unlist(lapply(Ns, function(i) sum(u.slice[i] < pi.prop))))
      Gs           <- seq_len(G)
      slice.ind    <- do.call(cbind, lapply(Gs, function(g, x=csi[g]) (u.slice < x)/x))
    
    # Cluster Labels
      psi          <- 1/psi.inv
      Sigma        <- lapply(Gs, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.res        <- sim.z(data=data, mu=mu[,Gs], Sigma=Sigma, Gseq=Gs, N=N, pi.prop=pi.prop[,Gs, drop=F], slice.ind=slice.ind, Q0=Q0s[Gs])
      z            <- z.res$z
      nn           <- tabulate(z, nbins=trunc.G)
      nn0          <- nn > 0
      z.ind        <- lapply(Gs, function(g) z == g)
      dat.g        <- lapply(Gs, function(g) data[z.ind[[g]],, drop=F])
      
    # Scores & Loadings
      c.data       <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0)   {
        f.tmp      <- lapply(Gs, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else matrix(, nr=0, nc=Q))
        FtF        <- lapply(Gs, function(g) if(nn0[g]) crossprod(f.tmp[[g]]))
        lmat[,,Gs] <- array(unlist(lapply(Gs, function(g) if(nn0[g]) matrix(unlist(lapply(Ps, function(j) sim.load(l.sigma=l.sigma, Q=Q, P=P, c.data=c.data[[g]][,j], f=f.tmp[[g]], Q1=Q1, 
                            psi.inv=psi.inv[,g][j], FtF=FtF[[g]], shrink=F)), use.names=F), nr=P, byrow=T) else sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F)), use.names=F), dim=c(P, Q, G))
        f          <- do.call(rbind, f.tmp)[obsnames,, drop=F]
      } else {
        f.tmp      <- lapply(Gs, function(g) f[z.ind[[g]],, drop=F])
      }
      
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.f        <- lapply(f.tmp, colSums)
      mu[,Gs]      <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], 
                              sum.f=sum.f[[g]], lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g]), mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero)))
                    
    # Uniquenesses
      psi.inv[,Gs] <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.psi.i(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, f=f.tmp[[g]], 
                              psi.beta=psi.beta, lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g])) else sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
      
    # Alpha
      if(MH.step)   {
        MH.alpha   <- sim.alpha(lower=MH.lower, upper=MH.upper, trunc.G=trunc.G, alpha=pi.alpha, Vs=Vs) 
        pi.alpha   <- MH.alpha$alpha
      }
      
      if(is.element(iter, iters))   {
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu 
        if(all(sw["f.sw"], Q0))    f.store[,,new.it]       <- f
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
        if(MH.step) {              rate[new.it]            <- MH.alpha$rate
                                   alpha.store[new.it]     <- pi.alpha }
                                   z.store[,new.it]        <- z 
                                   ll.store[new.it]        <- sum(z.res$log.likes)
                                   G.store[new.it]         <- sum(nn0)
      } 
    }
    returns        <- list(mu       = if(sw["mu.sw"])         mu.store,
                           f        = if(all(sw["f.sw"], Q0)) f.store, 
                           load     = if(all(sw["l.sw"], Q0)) load.store, 
                           psi      = if(sw["psi.sw"])        psi.store,
                           pi.prop  = if(sw["pi.sw"])         pi.store,
                           rate     = if(MH.step)             mean(rate),
                           alpha    = if(MH.step)             alpha.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store)
    return(returns)
  }