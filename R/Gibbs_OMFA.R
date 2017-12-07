###########################################################################
### Gibbs Sampler for Overfitted Mixtures of (Finite) Factor Analysers ####
###########################################################################

# Gibbs Sampler Function
  .gibbs_OMFA      <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu, mu,
                               sigma.l, burnin, thinning, uni.type, uni.prior,
                               sw, psi.alpha, psi.beta, verbose, cluster, ...) {

  # Define & initialise variables
    start.time     <- proc.time()
    sq_mat         <- if(P  > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix         <- base::matrix
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
    Q1             <- Q == 1
    uni            <- P == 1
    sw["s.sw"]     <- sw["s.sw"] && Q0
    sw["l.sw"]     <- sw["l.sw"] && Q0
    if(sw["mu.sw"])  {
      mu.store     <- array(0L,  dim=c(P, G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store    <- array(0L,  dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store   <- array(0L,  dim=c(P, Q, G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store    <- array(0L,  dim=c(P, G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store     <- matrix(0L, nrow=G, ncol=n.store)
    }
    z.store        <- matrix(0L, nrow=n.store, ncol=N)
    ll.store       <- vector("integer", n.store)
    err.z          <- zerr <- FALSE
    G.store        <- vector("integer", n.store)

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    nn.ind         <- which(nn > 0)
    pi.alpha       <- cluster$pi.alpha
    one.uni        <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv   <- switch(uni.type,   unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                         constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip    <- switch(uni.prior,  unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni)) {
      uni.shape    <- switch(uni.type,   constrained=N/2 + psi.alpha, single=(N * P)/2  + psi.alpha)
      V            <- switch(uni.type,   constrained=P, single=1)
    }
    psi.beta       <- switch(uni.prior,  isotropic=psi.beta[which.max(.ndeci(psi.beta))], psi.beta)
    pi.prop        <- c(cluster$pi.prop, vector("integer", G - length(cluster$pi.prop)))
    mu             <- cbind(mu,  vapply(seq_len(G - length(cluster$pi.prop)), function(g) .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P)))
    eta            <- .sim_eta_p(N=N, Q=Q)
    lmat           <- if(Q0) array(vapply(Gseq, function(g) .sim_load_p(Q=Q, P=P, sigma.l=sigma.l), numeric(P * Q)), dim=c(P, Q, G)) else array(0, dim=c(P, 0, G))
    psi.inv        <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    psi.inv        <- if(uni) t(psi.inv) else psi.inv
    if(Q0 &&   Q    < min(N - 1, Ledermann(P))) {
      for(g in which(nn     > P))   {
        fact       <- try(stats::factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]     <- fact$scores
          lmat[,,g]        <- unclass(fact$loadings)
          psi.inv[,g]      <- 1/fact$uniquenesses
        }
      }
    } else {
      psi.tmp      <- psi.inv
      if(isTRUE(one.uni)) {
        psi.inv[,] <- 1/switch(uni.type, constrained=Rfast::colVars(data), exp(mean(log(Rfast::colVars(data)))))
      } else   {
        psi.inv[,] <- 1/groupcolVars(data, ina=z)
      }
      inf.ind      <- is.infinite(psi.inv) | is.nan(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma        <- diag(1/sigma.l, Q)
    index          <- Order(nn, descending=TRUE)
    pi.prop        <- pi.prop[index]
    mu             <- mu[,index, drop=FALSE]
    lmat           <- lmat[,,index, drop=FALSE]
    psi.inv        <- psi.inv[,index, drop=FALSE]
    z              <- factor(z, labels=match(nn.ind, index))
    z              <- as.integer(levels(z))[z]
    G.non          <- G
    if(burnin       < 1)  {
      if(sw["mu.sw"])  mu.store[,,1]    <- mu
      if(sw["s.sw"])   eta.store[,,1]   <- eta
      if(sw["l.sw"])   load.store[,,,1] <- lmat
      if(sw["psi.sw"]) psi.store[,,1]   <- 1/psi.inv
      if(sw["pi.sw"])  pi.store[,1]     <- pi.prop
      z.store[1,]          <- z
      sigma                <- if(uni) lapply(Gseq, function(g) as.matrix(1/psi.inv[,g] + if(Q0) tcrossprod(as.matrix(lmat[,,g])) else 0)) else lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g]))
      log.probs            <- if(uni) vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log(pi.prop[g]), numeric(N)) else vapply(Gseq, function(g) { sigma <- if(Q0) sigma[[g]] else sq_mat(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q0) + log(pi.prop[g]) }, numeric(N))
      ll.store[1]          <- sum(rowLogSumExps(log.probs))
      G.store[1]           <- G.non
    }
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose   && iter  < burnin)  utils::setTxtProgressBar(pb, iter)
      storage      <- is.element(iter, iters)

    # Mixing Proportions & Re-ordering
      pi.prop      <- rDirichlet(G=G, alpha=pi.alpha, nn=nn)
      index        <- Order(nn, descending=TRUE)
      pi.prop      <- pi.prop[index]
      mu           <- mu[,index, drop=FALSE]
      lmat         <- lmat[,,index, drop=FALSE]
      psi.inv      <- psi.inv[,index, drop=FALSE]
      z            <- factor(z, labels=match(nn.ind, index))
      z            <- as.integer(levels(z))[z]

    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- if(uni) lapply(Gseq, function(g) as.matrix(psi[,g] + if(Q0) tcrossprod(as.matrix(lmat[,,g])) else 0)) else lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      log.pis      <- log(pi.prop)
      if(uni) {
        log.probs  <- vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log.pis[g], numeric(N))
      } else  {
        log.check  <- utils::capture.output(log.probs <- try(vapply(Gseq, function(g) dmvn(data, mu[,g], if(Q0) sigma[[g]] else sq_mat(sigma[[g]]), log=TRUE, isChol=!Q0) + log.pis[g], numeric(N)), silent=TRUE))
      }
      if(inherits(log.probs, "try-error")) {
       log.probs   <- vapply(Gseq, function(g) { sigma <- if(Q0) sigma[[g]] else sq_mat(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q0) + log.pis[g] }, numeric(N))
      }
      z            <- gumbel_max(probs=log.probs)
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      nn.ind       <- which(nn0)
      G.non        <- length(nn.ind)
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        eta.tmp    <- lapply(Gseq, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else .empty_mat(nc=Q))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) matrix(if(nn0[g]) vapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]],
                      Q1=Q1, EtE=EtE[[g]], psi.inv=psi.inv[,g][j]), numeric(Q)) else .sim_load_p(Q=Q, P=P, sigma.l=sigma.l), nrow=P, byrow=TRUE)), use.names=FALSE), dim=c(P, Q, G))
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else {
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat      <- lapply(Gseq, function(g) { S   <- c.data[[g]] - if(Q0) tcrossprod(eta.tmp[[g]], if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else 0; S * S } )
        psi.inv[,] <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,] <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], eta=eta.tmp[[g]], psi.beta=psi.beta,
                      P=P, Q0=Q0, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      }

    # Means
      sum.data     <- vapply(dat.g, colSums2, numeric(P))
      sum.data     <- if(uni) t(sum.data) else sum.data
      sum.eta      <- lapply(eta.tmp, colSums2)
      mu[,]        <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.zero=mu.zero, sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                             lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

      if(zerr && !err.z) {                                     warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(storage)  {
        if(verbose)   utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])             mu.store[,,new.it]      <- mu
        if(all(sw["s.sw"], Q0))     eta.store[,,new.it]     <- eta
        if(all(sw["l.sw"], Q0))     load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])            psi.store[,,new.it]     <- 1/psi.inv
        if(sw["pi.sw"])             pi.store[,new.it]       <- pi.prop
                                    z.store[new.it,]        <- as.integer(z)
                                    ll.store[new.it]        <- sum(rowLogSumExps(log.probs))
                                    G.store[new.it]         <- as.integer(G.non)
      }
    }
    if(verbose)       close(pb)
    Gmax           <- seq_len(max(as.integer(z.store)))
    mu.store       <- if(sw["mu.sw"])  tryCatch(mu.store[,Gmax,, drop=FALSE],    error=function(e) mu.store)
    load.store     <- if(sw["l.sw"])   tryCatch(load.store[,,Gmax,, drop=FALSE], error=function(e) load.store)
    psi.store      <- if(sw["psi.sw"]) tryCatch(psi.store[,Gmax,, drop=FALSE],   error=function(e) psi.store)
    pi.store       <- if(sw["pi.sw"])  tryCatch(pi.store[Gmax,, drop=FALSE],     error=function(e) pi.store)
    returns        <- list(mu        = if(sw["mu.sw"])         tryCatch(provideDimnames(mu.store,   base=list(varnames, "", ""),     unique=FALSE), error=function(e) mu.store),
                           eta       = if(all(sw["s.sw"], Q0)) tryCatch(provideDimnames(eta.store,  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                           load      = if(all(sw["l.sw"], Q0)) tryCatch(provideDimnames(load.store, base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                           psi       = if(sw["psi.sw"])        tryCatch(provideDimnames(psi.store,  base=list(varnames, "", ""),     unique=FALSE), error=function(e) psi.store),
                           pi.prop   = if(sw["pi.sw"])         pi.store,
                           z.store   = z.store,
                           ll.store  = ll.store,
                           G.store   = G.store,
                           time      = init.time)
    return(returns)
  }
