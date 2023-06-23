###########################################################
### Gibbs Sampler for Bayesian Infinite Factor Analysis ###
###########################################################

# Gibbs Sampler Function
  .gibbs_IFA     <- function(Q, data, iters, N, P, sigma.mu, mu, prop, truncated, uni.type,
                             uni.prior, psi.alpha, psi.beta, burnin, thinning, verbose,
                             sw, epsilon, mu.zero, nu1, nu2, adapt, active.crit, start.AGS, stop.AGS,
                             b0, b1, alpha.d1, alpha.d2, beta.d1, beta.d2, scaling, col.mean, ...) {

  # Define & initialise variables
    start.time   <- proc.time()
    matrix       <- base::matrix
    total        <- max(iters)
    if(verbose)     pb     <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    AGS.burn     <- total/5L
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    uni          <- P == 1
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
    post.mu      <- integer(P)
    post.psi     <- post.mu
    ll.store     <-
    Q.store      <- integer(n.store)
    Q.star       <- Q
    Q0           <- Q   > 0
    Q.large      <- Q.big  <- FALSE
    nu1.5        <- nu1 + 0.5
    P.5          <- P/2

    if(update.mu <- sw["u.sw"]) {
      mu.sigma   <- 1/sigma.mu
      mu.zero    <- as.numeric(mu.zero)
      mu.prior   <- mu.sigma * mu.zero
    } else mu[]  <- 0L
    uni.type     <- switch(EXPR=uni.type,  unconstrained=,               constrained="constrained", "single")
    .sim_psi_inv <- switch(EXPR=uni.type,  constrained=.sim_psi_u1,      single=.sim_psi_c1)
    .sim_psi_ip  <- switch(EXPR=uni.prior, unconstrained=.sim_psi_ipu,   isotropic=.sim_psi_ipc)
    psi.beta     <- switch(EXPR=uni.prior, isotropic=psi.beta[which.max(.ndeci(psi.beta))], psi.beta)
    uni.shape    <- switch(EXPR=uni.type,  constrained=N/2 + psi.alpha,  single=(N * P)/2 + psi.alpha)
    V            <- switch(EXPR=uni.type,  constrained=P,                single=1L)
    eta          <- .sim_eta_p(N=N, Q=Q)
    phi          <- .sim_phi_p(Q=Q, P=P, nu1=nu1, nu2=nu2)
    if(isTRUE(truncated)) {
     .sim_deltak <- .sim_deltaKT
     .rdelta     <- rltrgamma
     delta       <- c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_deltaPT(Q=Q, alpha=alpha.d2, beta=beta.d2))
    } else        {
      .rdelta    <- stats::rgamma
      delta      <- c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2))
    }
    tau          <- cumprod(delta)
    lmat         <- matrix(vapply(Pseq, function(j) .sim_load_ps(Q=Q, phi=phi[j,], tau=tau), numeric(Q)), nrow=P, byrow=TRUE)
    psi.inv      <- .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    psi.inv[]    <- 1/switch(EXPR=uni.type, constrained=colVars(data, center=col.mean, refine=FALSE, useNames=FALSE), max(colVars(data, center=col.mean, refine=FALSE, useNames=FALSE)))
    max.p        <- (psi.alpha  - 1)/psi.beta
    inf.ind      <- psi.inv > max(max.p)
    psi.inv[inf.ind]       <- switch(EXPR=uni.type, constrained=max.p, rep(max.p, P))[inf.ind]
    rm(max.p, inf.ind)
    sum.data     <- mu * N
    init.time    <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)) {
      if(verbose && iter    < burnin) utils::setTxtProgressBar(pb, iter)
      storage    <- is.element(iter,  iters)

    # Adaptation
      if(adapt   && all(iter >= start.AGS, iter < stop.AGS))      {
        if(stats::runif(1) < ifelse(iter < AGS.burn, 0.5, exp(-b0 - b1 * (iter - start.AGS)))) {
          switch(EXPR=active.crit, SC={
            if(Q0)   {
              SC <- .SC_crit(data, eta, lmat, prop)
              nonred      <- SC$nonred
              numred      <- SC$numred
            } else  {
              colvec      <- stats::runif(1)    <= prop
              numred      <- sum(colvec)
            }
          },        {
            colvec        <- if(Q0) (colSums2(abs(lmat) < epsilon, useNames=FALSE) / P) >= prop else stats::runif(1) <= prop
            numred        <- sum(colvec)
          })
          if(numred == 0)  {
            Q    <- Q + 1L
            Q.big   <- Q   > Q.star
            if(Q.big) {
              Q     <- Q.star
            } else {
              phi   <- cbind(phi,  .rgamma0(n=P, shape=nu1,      rate=nu2))
              delta <- c(delta,    .rdelta(n=1,  shape=alpha.d2, rate=beta.d2))
              tau   <- cumprod(delta)
              lmat  <- cbind(lmat, stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[,Q] * tau[Q])))
            }
          } else if(Q0)    {
            nonred  <- switch(EXPR=active.crit, BD=colvec == 0, nonred)
            Q       <- max(0L, Q - numred)
            phi     <- phi[,nonred, drop=FALSE]
            delta   <- delta[nonred]
            tau     <- cumprod(delta)
            lmat    <- lmat[,nonred, drop=FALSE]
          }
        }
      }
      Q0         <- Q  > 0
      Q1         <- Q == 1

    # Scores & Loadings
      c.data     <- sweep(data, 2L, mu, FUN="-", check.margin=FALSE)
      if(Q0) {
        eta      <- .sim_score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(vapply(Pseq, function(j) .sim_load_s(Q=Q, tau=tau, eta=eta, c.data=c.data[,j], Q1=Q1,
                           phi=phi[j,], psi.inv=psi.inv[j], EtE=crossprod(eta)), numeric(Q)), nrow=P, byrow=TRUE)
      } else {
        eta      <- .empty_mat(nr=N)
        lmat     <- .empty_mat(nr=P)
      }

    # Uniquenesses
      S.mat      <- c.data  - tcrossprod(eta, lmat)
      psi.inv[]  <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)

    # Means
      if(update.mu) {
        mu[]     <- .sim_mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums2(eta, useNames=FALSE), lmat=lmat, mu.prior=mu.prior)
      }

    # Shrinkage
      if(Q0) {
        load.2   <- lmat^2
        phi      <- .sim_phi(Q=Q, P=P, nu1.5=nu1.5, nu2=nu2, tau=tau, load.2=load.2)
        sum.term <- colSums2(phi * load.2, useNames=FALSE)
        for(k in seq_len(Q)) {
          delta[k]  <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[k], Q=Q, P.5=P.5, k=k,
                       tau.kq=tau[k:Q], sum.term.kq=sum.term[k:Q]) else .sim_delta1(Q=Q, P.5=P.5, tau=tau, sum.term=sum.term,
                       alpha.d1=ifelse(Q1, alpha.d2, alpha.d1), beta.d1=ifelse(Q1, beta.d2, beta.d1), delta.1=delta[1L])
          tau       <- cumprod(delta)
        }
      }

      if(Q.big && !Q.large && iter > burnin) {       cat("\n"); warning(paste0("\nQ has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star, "\n"), call.=FALSE)
        Q.large  <- TRUE
      }
      if(storage) {
        if(verbose) utils::setTxtProgressBar(pb, iter)
        new.it   <- which(iters == iter)
        psi      <- 1/psi.inv
        post.mu  <- post.mu  + mu/n.store
        post.psi <- post.psi + psi/n.store
        if(sw["mu.sw"])                          mu.store[,new.it] <- mu
        if(all(sw["s.sw"], Q0))      eta.store[,seq_len(Q),new.it] <- eta
        if(all(sw["l.sw"], Q0))     load.store[,seq_len(Q),new.it] <- lmat
        if(sw["psi.sw"])                        psi.store[,new.it] <- psi
                                                   Q.store[new.it] <- as.integer(Q)
                                                  ll.store[new.it] <- sum(dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + if(uni) psi else diag(psi), log=TRUE))
      }
    }
    if(verbose)     close(pb)
    Qmax         <- seq_len(max(Q.store))
    eta.store    <- if(sw["s.sw"])  tryCatch(eta.store[,Qmax,,  drop=FALSE], error=function(e) eta.store)
    load.store   <- if(sw["l.sw"])  tryCatch(load.store[,Qmax,, drop=FALSE], error=function(e) load.store)
    returns      <- list(mu       = if(sw["mu.sw"])  tryCatch(provideDimnames(mu.store,  base=list(varnames, ""), unique=FALSE), error=function(e) mu.store),
                         eta      = if(sw["s.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),        error=function(e) eta.store),  base=list(obsnames, "", ""), unique=FALSE), error=function(e) eta.store),
                         load     = if(sw["l.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),       error=function(e) load.store), base=list(varnames, "", ""), unique=FALSE), error=function(e) load.store),
                         psi      = if(sw["psi.sw"]) tryCatch(provideDimnames(psi.store, base=list(varnames, ""), unique=FALSE), error=function(e) psi.store),
                         post.mu  = tryCatch(stats::setNames(post.mu,  varnames),            error=function(e) post.mu),
                         post.psi = tryCatch(stats::setNames(post.psi, varnames),            error=function(e) post.psi),
                         ll.store = ll.store,
                         Q.store  = matrix(Q.store, nrow=1L),
                         time     = init.time)
    attr(returns, "Q.big") <- Q.large
      return(returns)
  }
