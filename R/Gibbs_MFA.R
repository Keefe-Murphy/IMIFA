#########################################################################
### Gibbs Sampler for (Finite) Mixtures of (Finite) Factor Analysers ####
#########################################################################

# Gibbs Sampler Function
  .gibbs_MFA       <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu, mu,
                               sigma.l, burnin, thinning, uni.type, uni.prior,
                               sw, psi.alpha, psi.beta, verbose, cluster, ...) {

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
    err.z          <- zerr <- FALSE
    ll.store       <- rep(0, n.store)

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    if(all(mu.zero  == 0)) {
      mu.zero      <- matrix(0, nrow=1, ncol=G)
      cluster$l.switch[1]  <- FALSE
    }
    if(length(mu.zero)  == 1) {
      mu.zero      <- matrix(mu.zero, nrow=1, ncol=G)
    }
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    pi.prop        <- cluster$pi.prop
    pi.alpha       <- cluster$pi.alpha
    .sim_psi_inv   <- switch(uni.type,  unconstrained=.sim_psi_iu,  isotropic=.sim_psi_ii)
    .sim_psi_ip    <- switch(uni.type,  unconstrained=.sim_psi_ipu, isotropic=.sim_psi_ipi)
    psi.beta       <- switch(uni.prior, isotropic=unique(round(psi.beta, min(nchar(psi.beta)))), psi.beta)
    if(length(psi.beta) == 1) {
      psi.beta     <- matrix(psi.beta, nrow=1, ncol=G)
    }
    mu0g           <- cluster$l.switch[1]
    psi0g          <- cluster$l.switch[2]
    label.switch   <- any(cluster$l.switch)
    eta            <- .sim_eta_p(N=N, Q=Q)
    lmat           <- lapply(Gseq, function(g) .sim_load_p(Q=Q, P=P, sigma.l=sigma.l))
    psi.inv        <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    if(Q0 && Q  < .ledermann(N, P)) {
      fact.ind     <- nn   <= P
      for(g in which(!fact.ind))    {
        fact       <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]     <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        }
      }
    } else     {
      psi.tmp      <- psi.inv
      psi.inv      <- vapply(Gseq, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma        <- diag(1/sigma.l, Q)
    lmat           <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      eta.store[,,1]       <- eta
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      sigma                <- lapply(Gseq, function(g) make.positive.definite(tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])))
      log.probs            <- vapply(Gseq, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      ll.store[1]          <- sum(gumbel_max(probs=log.probs, log.like=TRUE)$log.like)
    }
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose   && iter  < burnin) setTxtProgressBar(pb, iter)

    # Mixing Proportions
      pi.prop      <- rDirichlet(alpha=pi.alpha, G=G, nn=nn)

    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      log.check    <- capture.output(log.probs <- try(vapply(Gseq, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N)), silent=TRUE))
      if(inherits(log.probs, "try-error")) {
        log.probs  <- vapply(Gseq, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) make.positive.definite(sigma[[g]]) else make.positive.definite(sqrt(sigma[[g]])), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      }
      z.res        <- gumbel_max(probs=log.probs, log.like=TRUE)
      z            <- z.res$z
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        eta.tmp    <- lapply(Gseq, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else base::matrix(0, nrow=0, ncol=Q))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) if(nn0[g]) matrix(unlist(lapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]],
                      Q1=Q1, EtE=EtE[[g]], psi.inv=psi.inv[,g][j])), use.names=FALSE), nrow=P, byrow=TRUE) else .sim_load_p(Q=Q, P=P, sigma.l=sigma.l)), use.names=FALSE), dim=c(P, Q, G))
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else {
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
      }

    # Uniquenesses
      psi.inv      <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], eta=eta.tmp[[g]], psi.beta=psi.beta[,g],
                             P=P, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))

    # Means
      sum.data     <- vapply(dat.g, colSums, numeric(P))
      sum.eta      <- lapply(eta.tmp, colSums)
      mu           <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.zero=mu.zero[,g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g],  sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                             lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero[,g]), numeric(P))

    # Label Switching
      if(label.switch) {
        switch.lab <- .lab_switch(z.new=z, z.old=z.temp, Gs=Gseq)
        z          <- switch.lab$z
        z.perm     <- switch.lab$z.perm
        if(!identical(as.integer(z.perm), Gseq)) {
          mu       <- mu[,z.perm, drop=FALSE]
          lmat     <- lmat[,,z.perm, drop=FALSE]
          psi.inv  <- psi.inv[,z.perm, drop=FALSE]
          pi.prop  <- pi.prop[z.perm]
          nn       <- nn[z.perm]
         if(mu0g)  {
          mu.zero  <- mu.zero[,z.perm, drop=FALSE]
         }
         if(psi0g) {
          psi.beta <- psi.beta[,z.perm, drop=FALSE]
         }
        }
      }

      if(zerr && !err.z) {                                    warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(is.element(iter, iters))  {
        if(verbose)   setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu
        if(all(sw["s.sw"], Q0))    eta.store[,,new.it]     <- eta
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
                                   z.store[,new.it]        <- z
                                   ll.store[new.it]        <- sum(z.res$log.like)
      }
    }
    if(verbose)       close(pb)
    returns        <- list(mu       = if(sw["mu.sw"])         tryCatch(provideDimnames(mu.store,   base=list(varnames, "", ""),     unique=FALSE), error=function(e) mu.store),
                           eta      = if(all(sw["s.sw"], Q0)) tryCatch(provideDimnames(eta.store,  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                           load     = if(all(sw["l.sw"], Q0)) tryCatch(provideDimnames(load.store, base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                           psi      = if(sw["psi.sw"])        tryCatch(provideDimnames(psi.store,  base=list(varnames, "", ""),     unique=FALSE), error=function(e) psi.store),
                           pi.prop  = if(sw["pi.sw"])         pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           time     = init.time)
    attr(returns, "K")  <- mixFac_free(Q=Q, P=P, G=G, uni=uni.type)
    return(returns)
  }
