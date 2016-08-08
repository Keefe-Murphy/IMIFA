#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  gibbs.IMFA       <- function(Q, data, iters, N, P, G, mu.zero, pp, sigma.l, MH.step,
                               sigma.mu, burnin, thinning, mu, trunc.G, gen.slice, MH.lower,
                               psi.alpha, psi.beta, verbose, sw, cluster, MH.upper, ...) {
        
  # Define & initialise variables
    start.time     <- proc.time()
    total          <- max(iters)
    if(verbose)       pb   <- txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gs             <- seq_len(G)
    Ts             <- seq_len(trunc.G)
    Ps             <- seq_len(P)
    Ns             <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
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
    pi.prop        <- c(cluster$pi.prop, sim.pi(pi.alpha=pi.alpha, nn=rep(0, trunc.G), N=N, inf.G=TRUE, len=trunc.G)$pi.prop[-Gs])
    nn             <- tabulate(z, nbins=trunc.G)
    mu             <- cbind(mu, vapply(seq_len(trunc.G - G), function(g) sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P)))
    f              <- sim.f.p(N=N, Q=Q)
    lmat           <- lapply(Ts, function(t) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE))
    psi.inv        <- vapply(Ts, function(t) sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
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
      psi.inv      <- vapply(Ts, function(t) if(nn[t] > 1) 1/apply(data[z == t,, drop=FALSE], 2, var) else psi.tmp[,t], numeric(P))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[is.infinite(psi.inv)]
    }
    l.sigma        <- l.sigma * diag(Q)
    lmat           <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
    index          <- order(pi.prop, decreasing=TRUE)
    pi.prop        <- pi.prop[index]
    mu             <- mu[,index, drop=FALSE]
    lmat           <- lmat[,,index, drop=FALSE]
    psi.inv        <- psi.inv[,index, drop=FALSE]
    csi            <- pp * (1 - pp)^(Ts - 1)
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop/sum(pi.prop)
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu[,Gs], Gseq=Gs, N=N, pi.prop=pi.store[Gs,1], sigma=lapply(Gs,
                                  function(g) tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])), Q0=Q0s)$log.likes)
      G.store[1]           <- sum(nn > 0)
      if(MH.step)  {
        alpha.store[1]     <- pi.alpha 
      }
    }
    init.time      <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose   && iter  < burnin) setTxtProgressBar(pb, iter)
    
    # Mixing Proportions
      weights      <- sim.pi(pi.alpha=pi.alpha, nn=nn, N=N, inf.G=TRUE, len=trunc.G)
      pi.prop      <- weights$pi.prop
      if(MH.step)     Vs   <- weights$Vs
      
    # Slice Sampler
      if(!gen.slice) {
        index      <- order(pi.prop, decreasing=TRUE)
        pi.prop    <- csi  <- pi.prop[index]
        mu         <- mu[,index, drop=FALSE]
        lmat       <- lmat[,,index, drop=FALSE]
        psi.inv    <- psi.inv[,index, drop=FALSE]
      }
      u.slice      <- runif(N, 0, csi[z])
      G            <- max(vapply(Ns, function(i) sum(u.slice[i] < pi.prop), numeric(1)))
      Gs           <- seq_len(G)
      slice.ind    <- vapply(Gs, function(g, x=csi[g]) (u.slice < x)/x, numeric(N))
    
    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gs, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.res        <- sim.z(data=data, mu=mu[,Gs, drop=FALSE], sigma=sigma, Gseq=Gs, N=N, pi.prop=pi.prop[Gs], slice.ind=slice.ind, Q0=Q0s[Gs])
      z            <- z.res$z
      nn           <- tabulate(z, nbins=trunc.G)
      nn0          <- nn > 0
      z.ind        <- lapply(Gs, function(g) z == g)
      dat.g        <- lapply(Gs, function(g) data[z.ind[[g]],, drop=FALSE])
      
    # Scores & Loadings
      c.data       <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0)   {
        f.tmp      <- lapply(Gs, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else matrix(, nr=0, nc=Q))
        FtF        <- lapply(Gs, function(g) if(nn0[g]) crossprod(f.tmp[[g]]))
        lmat[,,Gs] <- array(unlist(lapply(Gs, function(g) if(nn0[g]) matrix(unlist(lapply(Ps, function(j) sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], f=f.tmp[[g]], Q1=Q1, FtF=FtF[[g]], 
                            P=P, psi.inv=psi.inv[,g][j], shrink=FALSE)), use.names=FALSE), nr=P, byrow=TRUE) else sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE)), use.names=FALSE), dim=c(P, Q, G))
        f          <- do.call(rbind, f.tmp)[obsnames,, drop=FALSE]
      } else {
        f.tmp      <- lapply(Gs, function(g) f[z.ind[[g]],, drop=FALSE])
      }
      
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.f        <- lapply(f.tmp, colSums)
      mu[,Gs]      <- vapply(Gs, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], sum.f=sum.f[[g]],
                             lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g]), mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P))
                    
    # Uniquenesses
      psi.inv[,Gs] <- vapply(Gs, function(g) if(nn0[g]) sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, f=f.tmp[[g]], psi.beta=psi.beta,
                             lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g])) else sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      
    # Alpha
      if(MH.step)   {
        MH.alpha   <- sim.alpha(lower=MH.lower, upper=MH.upper, trunc.G=trunc.G, alpha=pi.alpha, Vs=Vs) 
        pi.alpha   <- MH.alpha$alpha
      }
      
      if(is.element(iter, iters))   {
        if(verbose)   setTxtProgressBar(pb, iter)
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
    close(pb)
    Gmax           <- seq_len(max(as.numeric(z.store)))
    returns        <- list(mu       = if(sw["mu.sw"])         mu.store[,Gmax,, drop=FALSE],
                           f        = if(all(sw["f.sw"], Q0)) f.store, 
                           load     = if(all(sw["l.sw"], Q0)) load.store[,,Gmax,, drop=FALSE], 
                           psi      = if(sw["psi.sw"])        psi.store[,Gmax,, drop=FALSE],
                           pi.prop  = if(sw["pi.sw"])         pi.store[Gmax,, drop=FALSE],
                           rate     = if(MH.step)             mean(rate),
                           alpha    = if(MH.step)             alpha.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store,
                           time     = init.time)
    return(returns)
  }