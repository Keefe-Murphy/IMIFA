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
    total          <- max(iters)
    if(verbose)       pb   <- txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    Nseq           <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Gseq)
    iternames      <- paste0("Iteration", seq_len(n.store))
    Q0             <- Q  > 0
    Q0s            <- rep(Q0, G)
    Q1             <- Q == 1
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["s.sw"])   {
      eta.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)  <- list(obsnames, if(Q0) facnames, iternames)
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
    err.z          <- zerr <- FALSE
    G.store        <- rep(0, n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    
    mu.sigma       <- 1/sigma.mu
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    nn.ind         <- which(nn > 0)
    pi.alpha       <- cluster$pi.alpha
    pi.prop        <- cluster$pi.prop
    eta            <- sim.eta.p(N=N, Q=Q)
    lmat           <- lapply(Gseq, function(g) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE))
    psi.inv        <- vapply(Gseq, function(g) sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q0 && Q  < P - sqrt(P + Q)) {
      for(g in which(nn     > P))  {
        fact       <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]     <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } 
      }
    } else {
      psi.tmp      <- psi.inv
      psi.inv      <- vapply(Gseq, function(g) if(nn[g] > 1) 1/apply(data[z == g,, drop=FALSE], 2, var) else psi.tmp[,g], numeric(P))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma        <- diag(1/sigma.l, Q)
    lmat           <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
    index          <- order(nn, decreasing=TRUE)
    pi.prop        <- pi.prop[index]
    mu             <- mu[,index, drop=FALSE]
    lmat           <- lmat[,,index, drop=FALSE]
    psi.inv        <- psi.inv[,index, drop=FALSE]
    z              <- factor(z, labels=match(nn.ind, index))
    z              <- as.numeric(levels(z))[z]
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      eta.store[,,1]       <- eta
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, sigma=lapply(Gseq, function(g) 
                                        make.positive.definite(tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g]))), Q0=Q0s)$log.likes)
      G.store[1]           <- G
    }
    init.time      <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose   && iter  < burnin) setTxtProgressBar(pb, iter)
      
    # Mixing Proportions & Re-ordering
      pi.prop[]    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
      index        <- order(nn, decreasing=TRUE)
      pi.prop      <- pi.prop[index]
      mu           <- mu[,index, drop=FALSE]
      lmat         <- lmat[,,index, drop=FALSE]
      psi.inv      <- psi.inv[,index, drop=FALSE]
      z            <- factor(z, labels=match(nn.ind, index))
      z            <- as.numeric(levels(z))[z]
      
    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.log        <- capture.output({ z.res <- try(sim.z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s), silent=TRUE) })
      zerr         <- inherits(z.res, "try-error")
      if(zerr) {
        sigma      <- lapply(sigma, make.positive.definite)
        z.res      <- sim.z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s)
      }
      z            <- z.res$z
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      nn.ind       <- which(nn0)
      z.ind        <- lapply(Gseq, function(g) Nseq[z == g])
      dat.g        <- lapply(Gseq, function(g) data[z.ind[[g]],, drop=FALSE])
      
    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        eta.tmp    <- lapply(Gseq, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else base::matrix(0, nr=0, nc=Q))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) if(nn0[g]) matrix(unlist(lapply(Pseq, function(j) sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]], Q1=Q1, EtE=EtE[[g]], 
                            P=P, psi.inv=psi.inv[,g][j], shrink=FALSE)), use.names=FALSE), nr=P, byrow=TRUE) else sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE)), use.names=FALSE), dim=c(P, Q, G))
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else {
        eta.tmp    <- lapply(Gseq, function(g) eta[z.ind[[g]],, drop=FALSE])
      }
    
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.eta      <- lapply(eta.tmp, colSums)
      mu           <- vapply(Gseq, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], sum.eta=sum.eta[[g]], 
                             lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P))
      
    # Uniquenesses
      psi.inv      <- vapply(Gseq, function(g) if(nn0[g]) sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], eta=eta.tmp[[g]], psi.beta=psi.beta, 
                             P=P, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    
      if(zerr && !err.z) {                                    warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(is.element(iter, iters))   {
        if(verbose)   setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu 
        if(all(sw["s.sw"], Q0))    eta.store[,,new.it]     <- eta
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
                                   z.store[,new.it]        <- z 
                                   ll.store[new.it]        <- sum(z.res$log.likes)
                                   G.store[new.it]         <- sum(nn0)
      }
    }
    close(pb)
    Gmax           <- seq_len(max(as.numeric(z.store)))
    returns        <- list(mu       = if(sw["mu.sw"])         mu.store[,Gmax,, drop=FALSE],
                           eta      = if(all(sw["s.sw"], Q0)) eta.store, 
                           load     = if(all(sw["l.sw"], Q0)) load.store[,,Gmax,, drop=FALSE], 
                           psi      = if(sw["psi.sw"])        psi.store[,Gmax,, drop=FALSE],
                           pi.prop  = if(sw["pi.sw"])         pi.store[Gmax,, drop=FALSE],
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store,
                           time     = init.time)
    return(returns)
  }