#####################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Overfitted Case) ####
#####################################################################
  
# Gibbs Sampler Function
  gibbs.OMFA       <- function(Q, data, iters, N, P, G, mu.zero,
                               sigma.mu, sigma.l, burnin, mu, 
                               thinning, psi.alpha, psi.beta, 
                               sw, verbose, cluster, ...) {
        
  # Define & initialise variables
    start.time     <- proc.time()
    n.iters        <- round(max(iters), -1)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Gseq)
    iternames      <- paste0("Iteration", seq_len(n.store))
    Q0             <- Q > 0
    Q0s            <- rep(Q0, G)
    Q1             <- Q > 1
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store      <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store    <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store     <- matrix(0, nr=G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store        <- matrix(0, nr=N, nc=n.store)
    ll.store       <- rep(0, n.store)
    G.store        <- rep(0, n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    
    mu.sigma       <- 1/sigma.mu
    l.sigma        <- 1/sigma.l 
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    pi.alpha       <- cluster$pi.alpha
    pi.prop        <- cluster$pi.prop
    f              <- sim.f.p(N=N, Q=Q)
    lmat           <- lapply(Gseq, function(g) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE))
    psi.inv        <- do.call(cbind, lapply(Gseq, function(g) sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    if(Q0) {
      for(g in which(nn > 2.5 * Q))      {
        fact       <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } 
      }
    } else {
      psi.tmp      <- psi.inv
      psi.inv      <- do.call(cbind, lapply(Gseq, function(g) if(nn[g] > 1) 1/apply(data[z == g,, drop=FALSE], 2, var) else psi.tmp[,g]))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma        <- l.sigma * diag(Q)
    lmat           <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, sigma=lapply(Gseq,
                                  function(g) tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])), Q0=Q0s)$log.likes)
      G.store[1]           <- G
    }
    init.time      <- proc.time() - start.time
    
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
      pi.prop      <- sim.pi(pi.alpha=pi.alpha, nn=nn)
      
    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.res        <- sim.z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s)
      z            <- z.res$z
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      z.ind        <- lapply(Gseq, function(g) z == g)
      dat.g        <- lapply(Gseq, function(g) data[z.ind[[g]],, drop=FALSE])
      
    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        f.tmp      <- lapply(Gseq, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else matrix(, nr=0, nc=Q))
        FtF        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(f.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) if(nn0[g]) matrix(unlist(lapply(Pseq, function(j) sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], f=f.tmp[[g]], Q1=Q1, FtF=FtF[[g]], 
                            P=P, psi.inv=psi.inv[,g][j], shrink=FALSE)), use.names=FALSE), nr=P, byrow=TRUE) else sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE)), use.names=FALSE), dim=c(P, Q, G))
        f          <- do.call(rbind, f.tmp)[obsnames,, drop=FALSE]
      } else {
        f.tmp      <- lapply(Gseq, function(g) f[z.ind[[g]],, drop=FALSE])
      }
    
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.f        <- lapply(f.tmp, colSums)
      mu           <- do.call(cbind, lapply(Gseq, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], 
                              sum.f=sum.f[[g]], lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g]), mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero)))
    
    # Uniquenesses
      psi.inv      <- do.call(cbind, lapply(Gseq, function(g) if(nn0[g]) sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, f=f.tmp[[g]],
                              psi.beta=psi.beta, lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g])) else sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    
      if(is.element(iter, iters))   {
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu 
        if(all(sw["f.sw"], Q0))    f.store[,,new.it]       <- f
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
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
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store,
                           time     = init.time)
    return(returns)
  }