#####################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Overfitted Case) ####
#####################################################################

# Gibbs Sampler Function
  .gibbs_OMFA      <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu,
                               mu, sigma.l, burnin, thinning, sw, uni.type,
                               psi.alpha, psi.beta, verbose, cluster, ...) {

  # Define & initialise variables
    start.time     <- proc.time()
    total          <- max(iters)
    if(verbose)       pb   <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    Nseq           <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    Q0             <- Q  > 0
    Q0s            <- rep(Q0, G)
    Q1             <- Q == 1
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store    <- array(0, dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store    <- array(0, dim=c(P, G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store     <- matrix(0, nrow=G, ncol=n.store)
    }
    z.store        <- matrix(0, nrow=N, ncol=n.store)
    ll.store       <- rep(0, n.store)
    err.z          <- zerr <- FALSE
    G.store        <- rep(0, n.store)

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    nn.ind         <- which(nn > 0)
    pi.prop        <- cluster$pi.prop
    pi.alpha       <- cluster$pi.alpha
    .sim_psi.inv   <- switch(uni.type, unconstrained=.sim_psi.iu,  isotropic=.sim_psi.ii)
    .sim_psi.ip    <- switch(uni.type, unconstrained=.sim_psi.ipu, isotropic=.sim_psi.ipi)
    psi.beta       <- unique(round(psi.beta, min(nchar(psi.beta))))
    eta            <- .sim_eta.p(N=N, Q=Q)
    lmat           <- lapply(Gseq, function(g) .sim_load.p(Q=Q, P=P, sigma.l=sigma.l))
    psi.inv        <- vapply(Gseq, function(g) .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q0 && Q  < .ledermann(N, P)) {
      for(g in which(nn     > P))   {
        fact       <- try(stats::factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]     <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        }
      }
    } else {
      psi.tmp      <- psi.inv
      psi.inv      <- vapply(Gseq, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
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
      ll.store[1]          <- .sim_z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, sigma=lapply(Gseq, function(g)
                              corpcor::make.positive.definite(tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g]))), Q0=Q0s)$log.like
      G.store[1]           <- G
    }
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose   && iter  < burnin) utils::setTxtProgressBar(pb, iter)

    # Mixing Proportions & Re-ordering
      pi.prop      <- if(G == 1) 1 else .sim_pi(pi.alpha=pi.alpha, nn=nn, G)
      index        <- order(nn, decreasing=TRUE)
      pi.prop      <- pi.prop[index]
      mu           <- mu[,index, drop=FALSE]
      lmat         <- lmat[,,index, drop=FALSE]
      psi.inv      <- psi.inv[,index, drop=FALSE]
      z            <- factor(z, labels=match(nn.ind, index))
      z            <- as.numeric(levels(z))[z]

    # Cluster Labels
      psi          <- 1/psi.inv
      if(G > 1)  {
        sigma      <- lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
        z.log      <- utils::capture.output({ z.res <- try(.sim_z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s), silent=TRUE) })
        zerr       <- inherits(z.res, "try-error")
        if(zerr) {
          sigma    <- lapply(sigma, corpcor::make.positive.definite)
          z.res    <- .sim_z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s)
        }
        z          <- z.res$z
      } else     {
        z          <- rep(1, N)
      }
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      nn.ind       <- which(nn0)
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        eta.tmp    <- lapply(Gseq, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else base::matrix(0, nrow=0, ncol=Q))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) if(nn0[g]) matrix(unlist(lapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]],
                      Q1=Q1, EtE=EtE[[g]], psi.inv=psi.inv[,g][j])), use.names=FALSE), nrow=P, byrow=TRUE) else .sim_load.p(Q=Q, P=P, sigma.l=sigma.l)), use.names=FALSE), dim=c(P, Q, G))
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else {
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
      }

    # Means
      sum.data     <- vapply(dat.g, colSums, numeric(P))
      sum.eta      <- lapply(eta.tmp, colSums)
      mu           <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.zero=mu.zero, sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                             lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu.p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Uniquenesses
      psi.inv      <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], eta=eta.tmp[[g]], psi.beta=psi.beta,
                             P=P, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))

      if(zerr && !err.z) {                                    warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(is.element(iter, iters))   {
        if(verbose)    utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])             mu.store[,,new.it]     <- mu
        if(all(sw["s.sw"], Q0))     eta.store[,,new.it]    <- eta
        if(all(sw["l.sw"], Q0))     load.store[,,,new.it]  <- lmat
        if(sw["psi.sw"])            psi.store[,,new.it]    <- psi
        if(sw["pi.sw"])             pi.store[,new.it]      <- pi.prop
                                    z.store[,new.it]       <- z
                                    ll.store[new.it]       <- z.res$log.like
                                    G.store[new.it]        <- sum(nn0)
      }
    }
    close(pb)
    Gmax           <- seq_len(max(as.numeric(z.store)))
    returns        <- list(mu       = if(sw["mu.sw"])         provideDimnames(mu.store[,Gmax,, drop=FALSE],    base=list(varnames, "", ""),     unique=FALSE),
                           eta      = if(all(sw["s.sw"], Q0)) provideDimnames(eta.store,                       base=list(obsnames, "", ""),     unique=FALSE),
                           load     = if(all(sw["l.sw"], Q0)) provideDimnames(load.store[,,Gmax,, drop=FALSE], base=list(varnames, "", "", ""), unique=FALSE),
                           psi      = if(sw["psi.sw"])        provideDimnames(psi.store[,Gmax,, drop=FALSE],   base=list(varnames, "", ""),     unique=FALSE),
                           pi.prop  = if(sw["pi.sw"])         pi.store[Gmax,, drop=FALSE],
                           z.store  = z.store,
                           ll.store = ll.store,
                           G.store  = G.store,
                           time     = init.time)
    return(returns)
  }
