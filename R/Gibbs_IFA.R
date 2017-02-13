###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################

# Gibbs Sampler Function
  .gibbs_IFA     <- function(Q, data, iters, N, P, sigma.mu, mu, prop, uni.type,
                             psi.alpha, psi.beta, burnin, thinning, verbose, sw,
                             epsilon, mu.zero, nu, adapt, adapt.at, b0, b1,
                             alpha.d1, alpha.d2, beta.d1, beta.d2, nuplus1, ...) {

  # Define & initialise variables
    start.time   <- proc.time()
    total        <- max(iters)
    if(verbose)     pb     <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
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
    post.psi     <- post.mu
    ll.store     <- rep(0, n.store)
    cov.emp      <- Rfast::cova(as.matrix(data))
    cov.est      <- matrix(0, nrow=P, ncol=P)
    Q.star       <- Q
    Q.store      <- ll.store
    Q.large      <- Q.big  <- FALSE

    mu.sigma     <- 1/sigma.mu
    .sim_psi.inv <- switch(uni.type, unconstrained=.sim_psi.iu,  isotropic=.sim_psi.ii)
    .sim_psi.ip  <- switch(uni.type, unconstrained=.sim_psi.ipu, isotropic=.sim_psi.ipi)
    psi.beta     <- unique(round(psi.beta, min(nchar(psi.beta))))
    eta          <- .sim_eta.p(Q=Q, N=N)
    phi          <- .sim_phi.p(Q=Q, P=P, nu=nu, plus1=nuplus1)
    delta        <- c(.sim_delta.p(alpha=alpha.d1, beta=beta.d1), .sim_delta.p(Q=Q, alpha=alpha.d2, beta=beta.d2))
    tau          <- cumprod(delta)
    lmat         <- matrix(unlist(lapply(Pseq, function(j) .sim_load.ps(Q=Q, phi=phi[j,], tau=tau)), use.names=FALSE), nrow=P, byrow=TRUE)
    psi.inv      <- .sim_psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    if(all(Q      < .ledermann(N, P))) {
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
    sum.data     <- mu * N
    if(burnin     < 1) {
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
      Q0         <- Q  > 0
      Q1         <- Q == 1

    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q0) {
        eta      <- .sim_score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(unlist(lapply(Pseq, function(j) .sim_load.s(Q=Q, tau=tau, eta=eta, c.data=c.data[,j], Q1=Q1,
                           phi=phi[j,], psi.inv=psi.inv[j], EtE=crossprod(eta))), use.names=FALSE), nrow=P, byrow=TRUE)
      } else {
        eta      <- base::matrix(0, nrow=N, ncol=0)
        lmat     <- base::matrix(0, nrow=P, ncol=0)
      }

    # Means
      mu[]       <- .sim_mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums(eta), lmat=lmat, mu.zero=mu.zero)

    # Uniquenesses
      psi.inv    <- .sim_psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, c.data=c.data, eta=eta, lmat=lmat)

    # Local Shrinkage
      load.2     <- lmat * lmat
      phi        <- .sim_phi(Q=Q, P=P, nu=nu, tau=tau, load.2=load.2, plus1=nuplus1)

    # Global Shrinkage
      sum.term   <- diag(crossprod(phi, load.2))
      for(k in seq_len(Q)) {
        delta[k] <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[k], Q=Q, P=P, k=k,
                    tau.kq=tau[k:Q], sum.term.kq=sum.term[k:Q]) else .sim_delta1(Q=Q, P=P, tau=tau, sum.term=sum.term,
                    alpha.d1=ifelse(Q1, alpha.d2, alpha.d1), beta.d1=ifelse(Q1, beta.d2, beta.d1), delta.1=delta[1])
        tau      <- cumprod(delta)
      }

    # Adaptation
      if(all(adapt, iter > adapt.at)) {
        if(stats::runif(1) < ifelse(iter < burnin, 0.5, 1/exp(b0 + b1 * (iter - adapt.at)))) {
          colvec <- (if(Q0) colSums(abs(lmat) < epsilon) / P else 0) >= prop
          numred <- sum(colvec)
          if(numred == 0)  { # simulate extra columns from priors
            Q    <- Q + 1
            Q.big   <- Q > Q.star
            if(Q.big) {
              Q     <- Q.star
            } else {
              eta   <- cbind(eta, stats::rnorm(N))
              phi   <- cbind(phi, stats::rgamma(n=P, shape=nu + nuplus1, rate=nu))
              delta <- c(delta, stats::rgamma(n=1, shape=alpha.d2, rate=beta.d2))
              tau   <- cumprod(delta)
              lmat  <- cbind(lmat, stats::rnorm(n=P, mean=0, sd=sqrt(1/(phi[,Q] * tau[Q]))))
            }
          } else          { # remove redundant columns
            nonred  <- which(colvec == 0)
            Q       <- Q - numred
            eta     <- eta[,nonred, drop=FALSE]
            phi     <- phi[,nonred, drop=FALSE]
            delta   <- delta[nonred]
            tau     <- cumprod(delta)
            lmat    <- lmat[,nonred, drop=FALSE]
          }
        }
      }

      if(Q.big && !Q.large && iter > burnin) {       warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large  <- TRUE
      }
      if(is.element(iter, iters))  {
        if(verbose) utils::setTxtProgressBar(pb, iter)
        new.it   <- which(iters == iter)
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + sigma/n.store
        if(sw["mu.sw"])             mu.store[,new.it]              <- mu
        if(all(sw["s.sw"], Q0))     eta.store[,seq_len(Q),new.it]  <- eta
        if(all(sw["l.sw"], Q0))     load.store[,seq_len(Q),new.it] <- lmat
        if(sw["psi.sw"])            psi.store[,new.it]             <- psi
                                    Q.store[new.it]                <- Q
                                    ll.store[new.it]               <- sum(mvnfast::dmvn(X=data, mu=mu, sigma=sigma, log=TRUE))
      }
    }
    close(pb)
    Qmax         <- seq_len(max(Q.store))
    eta.store    <- if(sw["s.sw"])  tryCatch(eta.store[,Qmax,, drop=FALSE],  error=function(e) eta.store)
    load.store   <- if(sw["l.sw"])  tryCatch(load.store[,Qmax,, drop=FALSE], error=function(e) load.store)
    returns      <- list(mu       = if(sw["mu.sw"])  provideDimnames(mu.store,  base=list(varnames, ""), unique=FALSE),
                         eta      = if(sw["s.sw"])   provideDimnames(tryCatch(slam::as.simple_sparse_array(eta.store),  error=function(e) eta.store),  base=list(obsnames, "", ""), unique=FALSE),
                         load     = if(sw["l.sw"])   provideDimnames(tryCatch(slam::as.simple_sparse_array(load.store), error=function(e) load.store), base=list(varnames, "", ""), unique=FALSE),
                         psi      = if(sw["psi.sw"]) provideDimnames(psi.store, base=list(varnames, ""), unique=FALSE),
                         post.mu  = stats::setNames(post.mu,  varnames),
                         post.psi = stats::setNames(post.psi, varnames),
                         cov.emp  = provideDimnames(cov.emp,  base=list(varnames)),
                         cov.est  = provideDimnames(cov.est,  base=list(varnames)),
                         ll.store = ll.store,
                         Q.store  = matrix(Q.store, nrow=1),
                         time     = init.time)
    attr(returns, "Q.big") <- Q.large
    return(returns)
  }
