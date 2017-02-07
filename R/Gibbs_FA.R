################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################

# Gibbs Sampler Function
  .gibbs_FA      <- function(Q, data, iters, N, P, sigma.mu, mu, burnin,
                             thinning, uni.type, psi.alpha, psi.beta,
                             mu.zero, verbose, sw, sigma.l, ...) {

  # Define & initialise variables
    start.time   <- proc.time()
    total        <- max(iters)
    if(verbose)     pb     <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor", seq_len(Q))
    Q0           <- Q  > 0
    Q1           <- Q == 1
    dimnames(data)         <- NULL
    if(sw["mu.sw"])  {
      mu.store   <- matrix(0, nrow=P, ncol=n.store)
    }
    if(sw["s.sw"])   {
      eta.store  <- array(0, dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store  <- matrix(0, nrow=P, ncol=n.store)
    }
    post.mu      <- rep(0, P)
    post.psi     <- rep(0, P)
    ll.store     <- rep(0, n.store)
    cov.emp      <- Rfast::cova(as.matrix(data))
    cov.est      <- matrix(0, nrow=P, ncol=P)

    mu.sigma     <- 1/sigma.mu
    .sim_psi.inv <- switch(uni.type, unconstrained=.sim_psi.iu,  isotropic=.sim_psi.ii)
    .sim_psi.ip  <- switch(uni.type, unconstrained=.sim_psi.ipu, isotropic=.sim_psi.ipi)
    psi.beta     <- unique(round(psi.beta, min(nchar(psi.beta))))
    eta          <- .sim_eta.p(Q=Q, N=N)
    lmat         <- .sim_load.p(Q=Q, P=P, sigma.l=sigma.l)
    psi.inv      <- .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    if(all(Q0, Q  < .ledermann(N, P))) {
      fact       <- try(stats::factanal(data, factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
      if(!inherits(fact, "try-error")) {
        eta      <- fact$scores
        lmat     <- fact$loadings
        psi.inv  <- 1/fact$uniquenesses
      }
    } else {
      psi.tmp    <- psi.inv
      psi.inv    <- 1/Rfast::colVars(data)
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
      ll.store[1]          <- sum(mvnfast::dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + diag(1/psi.inv), log=TRUE))
    }
    init.time    <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose && iter    < burnin) utils::setTxtProgressBar(pb, iter)

    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q0) {
        eta      <- .sim_score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(unlist(lapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, Q1=Q1, c.data=c.data[,j],
                           eta=eta, psi.inv=psi.inv[j], EtE=crossprod(eta))), use.names=FALSE), nrow=P, byrow=TRUE)
      }

    # Means
      mu[]       <- .sim_mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums(eta), lmat=lmat, mu.zero=mu.zero)

    # Uniquenesses
      psi.inv    <- .sim_psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, c.data=c.data, eta=eta, lmat=lmat)

      if(is.element(iter, iters)) {
        if(verbose) utils::setTxtProgressBar(pb, iter)
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
                                 ll.store[new.it]     <- sum(mvnfast::dmvn(X=data, mu=mu, sigma=sigma, log=TRUE))
      }
    }
    close(pb)
    if(sw["s.sw"])               dimnames(eta.store)  <- list(obsnames, if(Q0) facnames, NULL)
    if(sw["l.sw"])               dimnames(load.store) <- list(varnames, if(Q0) facnames, NULL)
    returns   <- list(mu       = if(sw["mu.sw"])         provideDimnames(mu.store,   base=list(varnames, "")),
                      eta      = if(all(sw["s.sw"], Q0)) eta.store,
                      load     = if(all(sw["l.sw"], Q0)) load.store,
                      psi      = if(sw["psi.sw"])        provideDimnames(psi.store,  base=list(varnames, "")),
                      post.mu  = stats::setNames(post.mu,  varnames),
                      post.psi = stats::setNames(post.psi, varnames),
                      cov.emp  = provideDimnames(cov.emp, base=list(varnames, varnames)),
                      cov.est  = provideDimnames(cov.est, base=list(varnames, varnames)),
                      ll.store = ll.store,
                      time     = init.time)
    attr(returns, "K")        <- .dim(Q, P)
    return(returns)
  }
