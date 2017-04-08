############################################################################
### Gibbs Sampler for Infinite DP Mixtures of (Finite) Factor Analysers ####
############################################################################

# Gibbs Sampler Function
  .gibbs_IMFA        <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, learn.alpha, discount,  mu,
                                 a.hyper, sigma.mu, burnin, thinning, d.hyper, learn.d, uni.type, uni.prior, trunc.G,
                                 ind.slice, psi.alpha, psi.beta, verbose, sw, cluster, IM.lab.sw, zeta, kappa, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    total            <- max(iters)
    if(verbose)         pb   <- txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    Gs               <- seq_len(G)
    Ts               <- seq_len(trunc.G)
    Ps               <- seq_len(P)
    Ns               <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
    Q0               <- Q  > 0
    Q0s              <- rep(Q0, trunc.G)
    Q1               <- Q == 1
    if(sw["mu.sw"])  {
      mu.store       <- array(0, dim=c(P, trunc.G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store      <- array(0, dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store     <- array(0, dim=c(P, Q, trunc.G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store      <- array(0, dim=c(P, trunc.G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store       <- matrix(0, nrow=trunc.G, ncol=n.store)
    }
    z.store          <- matrix(0L, nrow=N, ncol=n.store)
    ll.store         <- rep(0,  n.store)
    acc1             <- acc2 <- FALSE
    err.z            <- zerr <- FALSE
    G.store          <- rep(0L, n.store)
    act.store        <- G.store
    if(learn.alpha) {
      alpha.store    <- ll.store
      alpha.shape    <- a.hyper[1]
      alpha.rate     <- a.hyper[2]
    }
    if(learn.d)     {
      d.store        <- ll.store
      d.shape1       <- d.hyper[1]
      d.shape2       <- d.hyper[2]
      d.rates        <- rep(0L, total)
      d.unif         <- d.shape1 == 1 & d.shape2 == 1
    } else d.rates   <- 1
    MH.step          <- any(discount  > 0, learn.d)
    if(MH.step)     {
      a.rates        <- rep(0L, total)
    } else a.rates   <- 1
    if(IM.lab.sw)   {
      lab.rate       <- matrix(0L, nrow=2, ncol=total)
    }
    mu.sigma         <- 1/sigma.mu
    sig.mu.sqrt      <- sqrt(sigma.mu)
    z                <- cluster$z
    pi.alpha         <- cluster$pi.alpha
    .sim_psi_inv     <- switch(uni.type,  unconstrained=.sim_psi_iu,  isotropic=.sim_psi_ii)
    .sim_psi_ip      <- switch(uni.type,  unconstrained=.sim_psi_ipu, isotropic=.sim_psi_ipi)
    psi.beta         <- switch(uni.prior, isotropic=unique(round(psi.beta, min(nchar(psi.beta)))), psi.beta)
    pi.prop          <- cluster$pi.prop
    nn               <- tabulate(z, nbins=trunc.G)
    mu               <- cbind(mu, vapply(seq_len(trunc.G - G), function(g) .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P)))
    eta              <- .sim_eta_p(N=N, Q=Q)
    lmat             <- lapply(Ts, function(t) .sim_load_p(Q=Q, P=P, sigma.l=sigma.l))
    psi.inv          <- vapply(Ts, function(t) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q0 && Q   < .ledermann(N, P)) {
      for(g in which(nn       > P))  {
        fact         <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]       <- fact$scores
          lmat[[g]]          <- fact$loadings
          psi.inv[,g]        <- 1/fact$uniquenesses
        }
      }
    } else {
      psi.tmp        <- psi.inv
      psi.inv[,Gs]   <- vapply(Gs, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
      inf.ind        <- is.infinite(psi.inv)
      psi.inv[inf.ind]       <- psi.tmp[inf.ind]
    }
    l.sigma          <- diag(1/sigma.l, Q)
    lmat             <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
    index            <- order(pi.prop, decreasing=TRUE)
    pi.prop          <- pi.prop[index]
    mu               <- mu[,index, drop=FALSE]
    lmat             <- lmat[,,index, drop=FALSE]
    psi.inv          <- psi.inv[,index, drop=FALSE]
    nn               <- nn[index]
    nn0              <- nn > 0
    nn.ind           <- which(nn0)
    G.non            <- length(nn.ind)
    z                <- factor(z, labels=match(nn.ind, index))
    z                <- as.integer(levels(z))[z]
    if(ind.slice) {
      ksi            <- (1 - rho) * rho^(Ts - 1)
      log.ksi        <- log(ksi)
    }
    slice.logs       <- c(- Inf, 0)
    if(burnin         < 1)  {
      mu.store[,,1]          <- mu
      eta.store[,,1]         <- eta
      load.store[,,,1]       <- lmat
      psi.store[,,1]         <- 1/psi.inv
      pi.store[,1]           <- pi.prop
      z.store[,1]            <- z
      sigma                  <- lapply(Gs, function(g) make.positive.definite(tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])))
      log.probs              <- vapply(Gs, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      ll.store[1]            <- sum(gumbel_max(probs=log.probs, log.like=TRUE)$log.like)
      G.store[1]             <- G.non
      act.store[1]           <- G
      if(learn.alpha) {
        alpha.store[1]       <- pi.alpha
      }
      if(learn.d)     {
        d.store[1]           <- discount
      }
    }
    init.time        <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose     && iter  < burnin)  setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

    # Mixing Proportions
      weights        <- if(ind.slice) .sim_pi_inf(alpha=pi.alpha, nn=nn[Gs], N=N, len=G, lseq=Gs, discount=discount) else .sim_pi_inf(alpha=pi.alpha, nn=nn, N=N, len=trunc.G, lseq=Ts, discount=discount)
      pi.prop        <- weights$pi.prop
      Vs             <- weights$Vs

    # Re-ordering & Slice Sampler
      index          <- order(pi.prop[Gs], decreasing=TRUE)
      pi.prop[Gs]    <- pi.prop[index]
      iVg            <- 1/Vs[G]
      Vs[Gs]         <- Vs[index]
      mu[,Gs]        <- mu[,index, drop=FALSE]
      lmat[,,Gs]     <- lmat[,,index, drop=FALSE]
      psi.inv[,Gs]   <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.integer(levels(z))[z]
      if(!ind.slice) {
        ksi          <- pi.prop
        log.ksi      <- log(ksi)
      }
      u.slice        <- runif(N, 0, ksi[z])
      G.old          <- G
      G              <- max(vapply(Ns, function(i) sum(u.slice[i] < ksi), integer(1L)))
      Gs             <- seq_len(G)
      if(ind.slice   && G.old < G) {
        newVs        <- if(discount == 0) rbeta(G - G.old, 1, pi.alpha) else rbeta(G - G.old, 1 - discount, pi.alpha + seq(G.old, G, 1) * discount)
        Vs           <- c(Vs,      newVs)
        newPis       <- newVs *    cumprod(c(pi.prop[G.old] * (iVg - 1), 1 - newVs[G.old - G]))
        pi.prop      <- c(pi.prop, newPis)
      } else pi.prop <- pi.prop[Gs]
      log.slice.ind  <- vapply(Gs, function(g) slice.logs[1 + (u.slice < ksi[g])] - log.ksi[g], numeric(N))

    # Cluster Labels
      psi            <- 1/psi.inv
      if(G > 1)  {
        sigma        <- lapply(Gs, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
        log.check    <- capture.output(log.probs  <- try(vapply(Gs, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N)), silent=TRUE))
        if(inherits(log.probs, "try-error")) {
          log.probs  <- vapply(Gs, function(g, Q=Q0s[g]) dmvn(data, mu[,g], if(Q) make.positive.definite(sigma[[g]]) else make.positive.definite(sqrt(sigma[[g]])), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
        }
        z.res        <- gumbel_max(probs=log.probs + log.slice.ind, log.like=TRUE, slice=TRUE)
        z            <- z.res$z
      } else     {
        z            <- rep(1, N)
      }
      nn             <- tabulate(z, nbins=trunc.G)
      nn0            <- nn > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)
      dat.g          <- lapply(Gs, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data         <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0)   {
        eta.tmp      <- lapply(Gs, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else base::matrix(0, nrow=0, ncol=Q))
        EtE          <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[,,Gs]   <- array(unlist(lapply(Gs, function(g) if(nn0[g]) matrix(unlist(lapply(Ps, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]],
                        EtE=EtE[[g]], Q1=Q1, psi.inv=psi.inv[,g][j])), use.names=FALSE), nrow=P, byrow=TRUE) else .sim_load_p(Q=Q, P=P, sigma.l=sigma.l)), use.names=FALSE), dim=c(P, Q, G))
        eta          <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else   {
        eta.tmp      <- lapply(Gs, function(g) eta[z == g,, drop=FALSE])
      }

    # Uniquenesses
      psi.inv[,Gs]   <- vapply(Gs, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, eta=eta.tmp[[g]], psi.beta=psi.beta,
                               lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))

    # Means
      sum.data       <- vapply(dat.g, colSums, numeric(P))
      sum.eta        <- lapply(eta.tmp, colSums)
      mu[,Gs]        <- vapply(Gs, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.zero=mu.zero, sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                               lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Alpha
      if(learn.alpha)      {
        if(discount  != 0) {
          MH.alpha   <- .sim_alpha_m(alpha=pi.alpha, discount=discount, alpha.shape=alpha.shape, alpha.rate=alpha.rate, N=N, G=G.non, zeta=zeta)
          pi.alpha   <- MH.alpha$alpha
          a.rate     <- MH.alpha$rate
        } else {
          pi.alpha   <- .sim_alpha_g(alpha=pi.alpha, shape=alpha.shape, rate=alpha.rate, G=G.non, N=N)
          a.rate     <- 1
        }
      }

    # Discount
      if(learn.d) {
        MH.disc      <- .sim_disc_mh(discount=discount, alpha=pi.alpha, disc.shape1=d.shape1, disc.shape2=d.shape2, N=N, G=G.non, kappa=kappa, unif=d.unif, nn=nn[nn0])
        discount     <- MH.disc$disc
        d.rate       <- MH.disc$rate
      }

    # Label Switching
      if(IM.lab.sw)   {
        if(G.non > 1) {
          move1      <- .lab_move1(nn.ind=nn.ind, pi.prop=pi.prop, nn=nn)
          acc1       <- move1$rate1
          if(acc1)    {
            sw1      <- move1$sw
            sw1x     <- c(sw1[2], sw1[1])
            nn[sw1]  <- nn[sw1x]
            nn0[sw1] <- nn0[sw1x]
            nn.ind   <- which(nn0)
            Vs[sw1]  <- Vs[sw1x]
            mu[,sw1] <- mu[,sw1x, drop=FALSE]
            lmat[,,sw1]      <- lmat[,,sw1x,   drop=FALSE]
            psi.inv[,sw1]    <- psi.inv[,sw1x, drop=FALSE]
            pi.prop[sw1]     <- pi.prop[sw1x]
            zsw1     <- z == sw1[1]
            z[z == sw1[2]]   <- sw1[1]
            z[zsw1]  <- sw1[2]
          }
        } else  acc1 <- FALSE
        if(G     > 1) {
          move2      <- .lab_move2(G=G, Vs=Vs, nn=nn)
          acc2       <- move2$rate2
          if(acc2)    {
            sw2      <- move2$sw
            sw2x     <- c(sw2[2], sw2[1])
            nn[sw2]  <- nn[sw2x]
            nn0[sw2] <- nn0[sw2x]
            nn.ind   <- which(nn0)
            mu[,sw2] <- mu[,sw2x, drop=FALSE]
            lmat[,,sw2]      <- lmat[,,sw2x,   drop=FALSE]
            psi.inv[,sw2]    <- psi.inv[,sw2x, drop=FALSE]
            pi.prop[sw2]     <- pi.prop[sw2x]
            zsw2     <- z == sw2[1]
            z[z == sw2[2]]   <- sw2[1]
            z[zsw2]  <- sw2[2]
          }
        } else  acc2 <- FALSE
      }
      if(zerr && !err.z) {                                       warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z        <- TRUE
      }
      if(MH.step)       a.rates[iter]                         <- a.rate
      if(learn.d)       d.rates[iter]                         <- d.rate
      if(IM.lab.sw)     lab.rate[,iter]                       <- c(acc1, acc2)
      if(storage)    {
        if(verbose)     setTxtProgressBar(pb, iter)
        new.it       <- which(iters == iter)
        if(sw["mu.sw"])               mu.store[,,new.it]      <- mu
        if(all(sw["s.sw"], Q0))       eta.store[,,new.it]     <- eta
        if(all(sw["l.sw"], Q0))       load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])              psi.store[,,new.it]     <- 1/psi.inv
        if(sw["pi.sw"])               pi.store[Gs,new.it]     <- pi.prop
        if(learn.alpha)               alpha.store[new.it]     <- pi.alpha
        if(learn.d)                   d.store[new.it]         <- discount
                                      z.store[,new.it]        <- as.integer(z)
                                      ll.store[new.it]        <- if(G > 1) sum(z.res$log.like) else sum(dmvn(X=data, mu=mu[,nn.ind], sigma=tcrossprod(as.matrix(lmat[,,nn.ind])) + diag(psi.store[,nn.ind,new.it]), log=TRUE))
                                      G.store[new.it]         <- as.integer(G.non)
                                      act.store[new.it]       <- as.integer(G)
      }
    }
    if(verbose)         close(pb)
    Gmax             <- seq_len(max(as.integer(z.store)))
    mu.store         <- if(sw["mu.sw"])  tryCatch(mu.store[,Gmax,, drop=FALSE],    error=function(e) mu.store)
    load.store       <- if(sw["l.sw"])   tryCatch(load.store[,,Gmax,, drop=FALSE], error=function(e) load.store)
    psi.store        <- if(sw["psi.sw"]) tryCatch(psi.store[,Gmax,, drop=FALSE],   error=function(e) psi.store)
    pi.store         <- if(sw["pi.sw"])  tryCatch(pi.store[Gmax,, drop=FALSE],     error=function(e) pi.store)
    returns          <- list(mu        = if(sw["mu.sw"])         tryCatch(provideDimnames(mu.store,   base=list(varnames, "", ""),     unique=FALSE), error=function(e) mu.store),
                             eta       = if(all(sw["s.sw"], Q0)) tryCatch(provideDimnames(eta.store,  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                             load      = if(all(sw["l.sw"], Q0)) tryCatch(provideDimnames(load.store, base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                             psi       = if(sw["psi.sw"])        tryCatch(provideDimnames(psi.store,  base=list(varnames, "", ""),     unique=FALSE), error=function(e) psi.store),
                             pi.prop   = if(sw["pi.sw"])         pi.store,
                             alpha     = if(learn.alpha)         alpha.store,
                             discount  = if(learn.d) {           if(sum(d.store == 0)/n.store > 0.5) as.simple_triplet_matrix(d.store) else d.store },
                             a.rate    = if(MH.step)             mean(a.rates) else a.rates,
                             d.rate    = if(learn.d)             mean(d.rates) else d.rates,
                             lab.rate  = if(IM.lab.sw)           setNames(rowmeans(lab.rate), c("Move1", "Move2")),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             G.store   = G.store,
                             act.store = act.store,
                             time      = init.time)
    return(returns)
  }
