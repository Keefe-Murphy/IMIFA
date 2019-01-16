################################################################
### Gibbs Sampler for Bayesian (Finite) Factor Analysis ########
################################################################

# Gibbs Sampler Function
  .gibbs_FA      <- function(Q, data, iters, N, P, sigma.mu, mu, burnin,
                             thinning, uni.type, uni.prior, psi.alpha, psi.beta,
                             sw, mu.zero, verbose, sigma.l, scaling, ...) {

  # Define & initialise variables
    start.time   <- proc.time()
    matrix       <- base::matrix
    total        <- max(iters)
    if(verbose)     pb     <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    Q0           <- Q  > 0
    Q1           <- Q == 1
    uni          <- P == 1
    sw["s.sw"]   <- sw["s.sw"] && Q0
    sw["l.sw"]   <- sw["l.sw"] && Q0
    dimnames(data)         <- NULL
    if(sw["mu.sw"])  {
      mu.store   <- matrix(0L, nrow=P, ncol=n.store)
    }
    if(sw["s.sw"])   {
      eta.store  <- array(0L,  dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store <- array(0L,  dim=c(P, Q, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store  <- matrix(0L, nrow=P, ncol=n.store)
    }
    post.mu      <- vector("integer", P)
    post.psi     <- post.mu
    ll.store     <- vector("integer", n.store)

    mu.sigma     <- 1/sigma.mu
    mu.zero      <- as.numeric(mu.zero)
    uni.type     <- switch(EXPR=uni.type,  unconstrained=,               constrained="constrained", "single")
    .sim_psi_inv <- switch(EXPR=uni.type,  constrained=.sim_psi_u1,      single=.sim_psi_c1)
    .sim_psi_ip  <- switch(EXPR=uni.prior, unconstrained=.sim_psi_ipu,   isotropic=.sim_psi_ipc)
    psi.beta     <- switch(EXPR=uni.prior, isotropic=psi.beta[which.max(.ndeci(psi.beta))], psi.beta)
    uni.shape    <- switch(EXPR=uni.type,  constrained=N/2 + psi.alpha,  single=(N * P)/2 + psi.alpha)
    V            <- switch(EXPR=uni.type,  constrained=P,                single=1L)
    eta          <- .sim_eta_p(Q=Q, N=N)
    lmat         <- matrix(.sim_load_p(Q=Q, P=P, sigma.l=sigma.l), nrow=P, ncol=Q)
    psi.inv      <- .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    psi.inv[]    <- 1/switch(EXPR=uni.type, constrained=.col_vars(data, avg=mu), .geom_mean(.col_vars(data, avg=mu)))
    max.p        <- (psi.alpha - 1)/switch(EXPR=uni.type, unconstrained=, constrained=psi.beta,
                                           min(do.call(psi_hyper, c(list(shape=psi.alpha, dat=data, covar=stats::cov(data)), list(...)$pots))))
    inf.ind      <- psi.inv > max(max.p)
    psi.inv[inf.ind]       <- switch(EXPR=uni.type, constrained=max.p, rep(max.p, P))[inf.ind]
    rm(max.p, inf.ind)
    l.sigma      <- diag(1/sigma.l, Q)
    sum.data     <- mu * N
    init.time    <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)) {
      if(verbose && iter    < burnin) utils::setTxtProgressBar(pb, iter)
      storage    <- is.element(iter,  iters)

    # Scores & Loadings
      c.data     <- sweep(data, 2L, mu, FUN="-", check.margin=FALSE)
      if(Q0) {
        eta      <- .sim_score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(vapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[,j], Q1=Q1,
                           eta=eta, psi.inv=psi.inv[j], EtE=crossprod(eta)), numeric(Q)), nrow=P, byrow=TRUE)
      }

    # Uniquenesses
      S.mat      <- c.data  - tcrossprod(eta, lmat)
      psi.inv[]  <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)

    # Means
      mu[]       <- .sim_mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums2(eta), lmat=lmat, mu.zero=mu.zero)

      if(storage) {
        if(verbose) utils::setTxtProgressBar(pb, iter)
        new.it   <- which(iters == iter)
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        if(sw["mu.sw"])             mu.store[,new.it]   <- mu
        if(all(sw["s.sw"], Q0))   eta.store[,,new.it]   <- eta
        if(all(sw["l.sw"], Q0))  load.store[,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,new.it]   <- psi
                                     ll.store[new.it]   <- sum(dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + if(uni) psi else diag(psi), log=TRUE))
      }
    }
    if(verbose)  close(pb)
    returns   <- list(mu       = if(sw["mu.sw"])           tryCatch(provideDimnames(mu.store,    base=list(varnames, ""),     unique=FALSE), error=function(e) mu.store),
                      eta      = if(all(sw["s.sw"], Q0))   tryCatch(provideDimnames(eta.store,   base=list(obsnames, "", ""), unique=FALSE), error=function(e) eta.store),
                      load     = if(all(sw["l.sw"], Q0))   tryCatch(provideDimnames(load.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) load.store),
                      psi      = if(sw["psi.sw"])          tryCatch(provideDimnames(psi.store,   base=list(varnames, ""),     unique=FALSE), error=function(e) psi.store),
                      post.mu  = tryCatch(stats::setNames(post.mu,  varnames),            error=function(e) post.mu),
                      post.psi = tryCatch(stats::setNames(post.psi, varnames),            error=function(e) post.psi),
                      ll.store = ll.store,
                      time     = init.time)
    attr(returns, "K")        <- PGMM_dfree(Q=Q, P=P, method=switch(EXPR=uni.type, constrained="UCU", single="UCC"))
      return(returns)
  }
