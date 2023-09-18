############################################################################
### Gibbs Sampler for Infinite DP Mixtures of (Finite) Factor Analysers ####
############################################################################

# Gibbs Sampler Function
  .gibbs_IMFA        <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, learn.alpha, discount, mu, tune.zeta,
                                 a.hyper, sigma.mu, burnin, thinning, d.hyper, learn.d, uni.type, uni.prior, trunc.G, col.mean,
                                 ind.slice, psi.alpha, psi.beta, verbose, sw, cluster, IM.lab.sw, zeta, kappa, thresh, exchange, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    sq_mat           <- if(P  > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix           <- base::matrix
    total            <- max(iters)
    if(verbose)         pb   <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    Gs               <- seq_len(G)
    Ts               <- seq_len(trunc.G)
    Ps               <- seq_len(P)
    Ns               <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
    Q0               <- Q  > 0
    Q1               <- Q == 1
    uni              <- P == 1
    sw["s.sw"]       <- sw["s.sw"] && Q0
    sw["l.sw"]       <- sw["l.sw"] && Q0
    if(sw["mu.sw"])  {
      mu.store       <- array(0L,  dim=c(P, trunc.G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store      <- array(0L,  dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store     <- array(0L,  dim=c(P, Q, trunc.G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store      <- array(0L,  dim=c(P, trunc.G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store       <- matrix(0L, nrow=trunc.G, ncol=n.store)
    }
    z.store          <- matrix(0L, nrow=n.store, ncol=N)
    ll.store         <-
    G.store          <- integer(n.store)
    acc1             <- acc2 <- FALSE
    err.z            <- zerr <- FALSE
    act.store        <- G.store
    pi.alpha         <- cluster$pi.alpha
    if(learn.alpha) {
      alpha.store    <- ll.store
      alpha.shape    <- a.hyper[1L]
      alpha.rate     <- a.hyper[2L]
    }
    if(learn.d)     {
      d.store        <- ll.store
      d.shape1       <- d.hyper[1L]
      d.shape2       <- d.hyper[2L]
      d.rates        <- integer(total)
      d.unif         <- d.shape1 == 1   && d.shape2 == 1
      .sim_disc_mh   <- if(!learn.alpha && pi.alpha == 0) .sim_d_slab else .sim_d_spike
    } else d.rates   <- 1L
    Dneg             <- !learn.d        && discount  < 0
    MH.step          <- any(discount  > 0, learn.d) && learn.alpha
    if(MH.step)     {
      a.rates        <- integer(total)
    } else a.rates   <- 1L
    if(IM.lab.sw)   {
      lab.rate       <- matrix(0L, nrow=2L, ncol=total)
    }
    abs.disc         <- abs(discount)
    d.count          <- 0L
    avgzeta          <- zeta
    heat             <- tune.zeta$heat
    lambda           <- tune.zeta$lambda
    target           <- tune.zeta$target
    zeta.tune        <- tune.zeta$do
    startz           <- tune.zeta$start.zeta
    stopz            <- tune.zeta$stop.zeta
    mu.sigma         <- 1/sigma.mu
    mu.prior         <- mu.sigma * mu.zero
    sig.mu.sqrt      <- sqrt(sigma.mu)
    l.sigma          <- diag(1/sigma.l, Q)
    sig.l.sqrt       <- sqrt(sigma.l)
    z                <- cluster$z
    nn               <- tabulate(z, nbins=trunc.G)
    nn0              <- nn > 0
    nn.ind           <- which(nn > 0)
    G.non            <- length(nn.ind)
    one.uni          <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv     <- switch(EXPR=uni.type,  unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                               constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip      <- switch(EXPR=uni.prior, unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni)) {
      uni.shape      <- switch(EXPR=uni.type,  constrained=N/2 + psi.alpha, single=(N * P)/2 + psi.alpha)
      V              <- switch(EXPR=uni.type,  constrained=P, single=1L)
    }
    psi.beta         <- switch(EXPR=uni.prior, isotropic=psi.beta[which.max(.ndeci(psi.beta))], psi.beta)
    pi.prop          <- cluster$pi.prop
    mu.tmp           <- vapply(seq_len(trunc.G - G), function(g) .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))
    mu               <- cbind(mu, if(uni) t(mu.tmp) else mu.tmp)
    eta              <- .sim_eta_p(N=N, Q=Q)
    if(Q0)            {
      eta            <- .sim_eta_p(N=N, Q=Q)
      lmat           <- array(vapply(Ts, function(t) .sim_load_p(Q=Q, P=P, sig.l.sqrt=sig.l.sqrt), numeric(P * Q)), dim=c(P, Q, trunc.G))
    } else            {
      eta            <- .empty_mat(nr=N)
      lmat           <- array(0L, dim=c(P, 0, trunc.G))
    }
    if(isTRUE(one.uni)) {
      psi.inv        <- matrix(.sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), nrow=P, ncol=trunc.G)
    } else psi.inv   <- replicate(trunc.G, .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), simplify="array")
    psi.inv          <- if(uni) t(psi.inv) else psi.inv
    if(isTRUE(one.uni)) {
      psi.inv[]      <- 1/switch(EXPR=uni.type, constrained=colVars(data, center=col.mean, refine=FALSE, useNames=FALSE), max(colVars(data, center=col.mean, refine=FALSE, useNames=FALSE)))
    } else   {
      tmp.psi        <- (nn[nn0] - 1L)/pmax(rowsum(data^2, z) - rowsum(data, z)^2/nn[nn0], 0L)
      tmp.psi        <- switch(EXPR=uni.type, unconstrained=t(tmp.psi), matrix(Rfast::rowMaxs(tmp.psi, value=TRUE), nrow=P, ncol=G, byrow=TRUE))
      psi.inv[,nn     > 1]   <- tmp.psi[!is.nan(tmp.psi)]
      rm(tmp.psi)
    }
    max.p            <- (psi.alpha  - 1)/psi.beta
    inf.ind          <- psi.inv  > max(max.p)
    psi.inv[inf.ind] <- matrix(max.p, nrow=P, ncol=trunc.G)[inf.ind]
    rm(max.p, inf.ind)
    if(ind.slice) {
      ksi            <- (1 - rho)   * rho^(Ts  - 1L)
      log.ksi        <- log(ksi)
      slinf          <- rep(-Inf, N)
    } else slinf     <- c(-Inf,  0L)
    if((noLearn      <-
      (isFALSE(learn.alpha)     &&
       isFALSE(learn.d)))       &&
       isTRUE(thresh))           {
      TRX            <- .slice_threshold(N, pi.alpha, discount, MPFR=pi.alpha == 0)
    }
    init.time        <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total))  {
      if(verbose     && iter  < burnin)  utils::setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

    # Mixing Proportions & Re-ordering
      if(!exchange)   {
        Vs           <- .sim_vs_inf(alpha=pi.alpha, nn=nn[Gs], N=N, discount=discount, len=G, lseq=Gs)
        pi.prop      <- .sim_pi_inf(Vs, len=G)
        prev.prod    <- 1 - sum(pi.prop)
        prev.prod    <- ifelse(prev.prod < 0, pi.prop[G] * (1/Vs[G] - 1), prev.prod)
        index        <- order(pi.prop, decreasing=TRUE)
      } else          {
        piX          <- .sim_pi_infX(nn=nn[nn0], Kn=G.non, G=G, alpha=pi.alpha, discount=discount)
        pi.prop      <- piX$pi.prop
        prev.prod    <- piX$prev.prod
        index        <- piX$index
      }
      GI             <- which(Gs[index] == G)
      pi.prop        <- pi.prop[index]
      Vs             <- if(!exchange) Vs[index] else integer(G)
      mu[,Gs]        <- mu[,index, drop=FALSE]
      lmat[,,Gs]     <- lmat[,,index, drop=FALSE]
      psi.inv[,Gs]   <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.integer(levels(z))[z]
      nn[Gs]         <- nn[index]
      nn0[Gs]        <- nn0[index]

    # Scores & Loadings
      dat.g          <- lapply(Gs, function(g) data[z == g,, drop=FALSE])
      c.data         <- lapply(Gs, function(g) sweep(dat.g[[g]], 2L, mu[,g], FUN="-", check.margin=FALSE))
      if(Q0)     {
        eta.tmp      <- lapply(Gs, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else .empty_mat(nc=Q))
        EtE          <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[,,Gs]   <- array(unlist(lapply(Gs, function(g) matrix(if(nn0[g]) vapply(Ps, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]], Q1=Q1,
                        EtE=EtE[[g]], psi.inv=psi.inv[,g][j]), numeric(Q)) else .sim_load_p(Q=Q, P=P, sig.l.sqrt=sig.l.sqrt), nrow=P, byrow=TRUE)), use.names=FALSE), dim=c(P, Q, G))
        eta          <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else     {
        eta.tmp      <- lapply(Gs, function(g) eta[z == g,, drop=FALSE])
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat        <- lapply(Gs, function(g) { S   <- c.data[[g]] - if(Q0) tcrossprod(eta.tmp[[g]], if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else 0L; S^2 } )
        psi.inv[,]   <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,Gs] <- vapply(Gs, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, eta=eta.tmp[[g]], psi.beta=psi.beta,
                               Q0=Q0, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      }

    # Means
      sum.data       <- vapply(dat.g, colSums2, useNames=FALSE, numeric(P))
      sum.data       <- if(uni) t(sum.data) else sum.data
      sum.eta        <- lapply(eta.tmp, colSums2, useNames=FALSE)
      mu[,Gs]        <- vapply(Gs, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.prior=mu.prior, sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                               lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Slice Sampler
      if(thresh      &&
         !noLearn)    {
        TRX          <- .slice_threshold(N, pi.alpha, discount, MPFR=pi.alpha == 0)
      }
      if(!ind.slice)  {
        ksi          <- if(thresh) pmin(pi.prop, TRX) else pi.prop
        log.ksi      <- log(ksi)
      }
      u.slice        <- stats::runif(N, 0L, ksi[z])
      min.u          <- min(u.slice)
      G.old          <- G
      if(ind.slice)   {
        G.new        <- sum(min.u < ksi)
        G.trunc      <- G.new < G.old
        while(G  < G.new && (Vs[GI] != 1 || pi.prop[GI]  != 0)) {
          newVs      <- .sim_vs_inf(alpha=pi.alpha, discount=discount, len=1L, lseq=G + 1L)
          Vs         <- c(Vs,      newVs)
          pi.prop    <- c(pi.prop, newVs  * prev.prod)
          GI    <- G <- G + 1L
          prev.prod  <- prev.prod * (1    - newVs)
          prev.prod  <- ifelse(prev.prod  < 0, pi.prop[G] * (1/Vs[G] - 1), prev.prod)
        }
        G            <- ifelse(G.trunc, G.new, G)
        Gs           <- seq_len(G)
        if(G.trunc)   {
          pi.prop    <- pi.prop[Gs]
          Vs         <- Vs[Gs]
        }
      } else     {
        cum.pi       <- sum(pi.prop)
        u.max        <- 1 - min.u
        G.trunc      <- cum.pi > u.max
        while(cum.pi  < u.max && trunc.G > G && (pi.prop[GI] != 0 || ifelse(exchange, prev.prod > min.u, Vs[GI] != 1))) {
          newVs      <- .sim_vs_inf(alpha=pi.alpha, discount=discount, len=1L, lseq=G + 1L)
          Vs         <- c(Vs,      newVs)
          newPis     <- newVs  *   prev.prod
          pi.prop    <-
          ksi        <- c(pi.prop, newPis)
          log.ksi    <- c(log.ksi, log(newPis))
          cum.pi     <- cum.pi +   newPis
          GI    <- G <- G + 1L
          prev.prod  <- 1 - cum.pi
          prev.prod  <- ifelse(prev.prod < 0, pi.prop[G] * (1/Vs[G] - 1), prev.prod)
        }
        G            <- ifelse(G.trunc, which.max(cumsum(pi.prop) > u.max), G)
        Gs           <- seq_len(G)
        if(G.trunc)   {
          pi.prop    <- ksi   <-   pi.prop[Gs]
          Vs         <- Vs[Gs]
        }
      }
      if(thresh)      {
        ksi          <- pmax(pi.prop, TRX)
        log.ksi      <- log(ksi)
      }

    # Cluster Labels
      if(G > 1)  {
        psi          <- 1/psi.inv
        sigma        <- if(uni) lapply(Gs, function(g) as.matrix(psi[,g] + if(Q0) tcrossprod(as.matrix(lmat[,,g])) else 0L)) else lapply(Gs, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
        if(ind.slice) {
          log.pixi   <- log(pi.prop) - log.ksi[Gs]
          log.probs  <- vapply(Gs, function(g) { slinf[u.slice < ksi[g]] <- log.pixi[g]; slinf }, numeric(N))
        } else   {
          log.probs  <- vapply(Gs, function(g) slinf[1L + (u.slice < ksi[g])], numeric(N))
        }
        fin.probs    <- is.finite(log.probs)
        if(uni)  {
          log.probs  <- vapply(Gs, function(g, LP=log.probs[,g], FP=fin.probs[,g]) { LP[FP] <- stats::dnorm(data[FP,], mu[,g], sq_mat(sigma[[g]]), log=TRUE) + LP[FP]; LP }, numeric(N))
        } else   {
          probs.log  <- log.probs
          log.probs  <- try(vapply(Gs, function(g, LP=log.probs[,g], FP=fin.probs[,g]) { LP[FP] <- dmvn(data[FP,], mu[,g], sigma[[g]], log=TRUE,  isChol=!Q) + LP[FP]; LP }, numeric(N)), silent=TRUE)
          if(zerr    <- inherits(log.probs, "try-error")) {
           log.probs <- vapply(Gs, function(g, LP=probs.log[,g], FP=fin.probs[,g], Q=Q0[g]) { sigma <- if(Q) is.posi_def(sigma[[g]], make=TRUE)$X.new else sq_mat(sigma[[g]]); LP[FP] <- dmvn(data[FP,], mu[,g], sigma, log=TRUE, isChol=!Q) + LP[FP]; LP }, numeric(N))
          }
          rm(probs.log)
        }
        z            <- gumbel_max(probs=log.probs, slice=TRUE)
      } else     {
        z            <- rep(1L, N)
      }
      nn             <- tabulate(z, nbins=trunc.G)
      nn0            <- nn > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)

    # Alpha
      if(learn.alpha)      {
        if((non0d    <- discount != 0)  && Dneg) {
          pi.alpha   <- G  * abs.disc
        } else if(isTRUE(non0d))  {
          MH.alpha   <- .sim_alpha_m(alpha=pi.alpha, discount=discount, alpha.shape=alpha.shape, alpha.rate=alpha.rate, N=N, G=G.non, zeta=zeta)
          pi.alpha   <- MH.alpha$alpha
          a.rate     <- MH.alpha$rate
          if(isTRUE(zeta.tune)) {
            d.count  <- d.count + non0d
            if(iter  >= startz &&
               iter   < stopz)  {
              zeta   <- .tune_zeta(zeta=zeta, time=d.count, l.rate=MH.alpha$l.prob, heat=heat, target=target, lambda=lambda)
            }
            avgzeta  <- c(avgzeta, zeta)
          }
        } else {
          pi.alpha   <- .sim_alpha_g(alpha=pi.alpha, shape=alpha.shape, rate=alpha.rate, G=G.non, N=N)
          a.rate     <- 1L
        }
      }

    # Discount
      if(learn.d) {
        MH.disc      <- .sim_disc_mh(discount=discount, disc.shape1=d.shape1, disc.shape2=d.shape2, G=G.non, unif=d.unif, nn=nn[nn0], alpha=pi.alpha, kappa=kappa)
        discount     <- MH.disc$disc
        d.rate       <- MH.disc$rate
      }

    # Label Switching
      if(IM.lab.sw)   {
        if(G.non > 1) {
          move1      <- .lab_move1(nn.ind=nn.ind, pi.prop=pi.prop, nn=nn)
          if((acc1   <- move1$rate1)) {
            sw1      <- move1$sw
            sw1x     <- c(sw1[2L], sw1[1L])
            nn[sw1]  <- nn[sw1x]
            nn0[sw1] <- nn0[sw1x]
            nn.ind   <- which(nn0)
            Vs[sw1]  <- Vs[sw1x]
            mu[,sw1] <- mu[,sw1x, drop=FALSE]
            lmat[,,sw1]      <- lmat[,,sw1x,   drop=FALSE]
            psi.inv[,sw1]    <- psi.inv[,sw1x, drop=FALSE]
            pi.prop[sw1]     <- pi.prop[sw1x]
            zsw1     <- z == sw1[1L]
            z[z == sw1[2L]]  <- sw1[1L]
            z[zsw1]  <- sw1[2L]
          }
        } else  acc1 <- FALSE
        if(G     > 1) {
          move2      <- .lab_move2(G=G, Vs=Vs, nn=nn)
          if((acc2   <- move2$rate2)) {
            sw2      <- move2$sw
            sw2x     <- c(sw2[2L], sw2[1L])
            nn[sw2]  <- nn[sw2x]
            nn0[sw2] <- nn0[sw2x]
            nn.ind   <- which(nn0)
            mu[,sw2] <- mu[,sw2x, drop=FALSE]
            lmat[,,sw2]      <- lmat[,,sw2x,   drop=FALSE]
            psi.inv[,sw2]    <- psi.inv[,sw2x, drop=FALSE]
            pi.prop[sw2]     <- pi.prop[sw2x]
            zsw2     <- z == sw2[1L]
            z[z == sw2[2L]]  <- sw2[1L]
            z[zsw2]  <- sw2[2L]
          }
        } else  acc2 <- FALSE
      }
      if(zerr && !err.z) {                                       cat("\n"); warning("\nAlgorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices\n", call.=FALSE, immediate.=TRUE)
        err.z        <- TRUE
      }
      if(MH.step)       a.rates[iter]                         <- a.rate
      if(learn.d)       d.rates[iter]                         <- d.rate
      if(IM.lab.sw)     lab.rate[,iter]                       <- c(acc1, acc2)
      if(storage)    {
        if(verbose)     utils::setTxtProgressBar(pb, iter)
        new.it       <- which(iters == iter)
        if(sw["mu.sw"])                  mu.store[,,new.it]   <- mu
        if(all(sw["s.sw"], Q0))         eta.store[,,new.it]   <- eta
        if(all(sw["l.sw"], Q0))       load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])                psi.store[,,new.it]   <- 1/psi.inv
        if(sw["pi.sw"])                 pi.store[Gs,new.it]   <- pi.prop
        if(learn.alpha)                 alpha.store[new.it]   <- pi.alpha
        if(learn.d)                         d.store[new.it]   <- discount
                                           z.store[new.it,]   <- as.integer(z)
                                           ll.store[new.it]   <- if(G > 1) sum(rowLogSumExps(log.probs, useNames=FALSE)) else sum(dmvn(X=data, mu=mu[,nn.ind], sigma=tcrossprod(as.matrix(lmat[,,nn.ind])) + diag(psi.store[,nn.ind,new.it]), log=TRUE))
                                            G.store[new.it]   <- as.integer(G.non)
                                          act.store[new.it]   <- as.integer(G)
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
                             a.rate    = ifelse(MH.step,         mean(a.rates), a.rates),
                             d.rate    = ifelse(learn.d,         mean(d.rates), d.rates),
                             lab.rate  = if(IM.lab.sw)           stats::setNames(rowMeans2(lab.rate, useNames=FALSE), c("Move1", "Move2")),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             G.store   = G.store,
                             act.store = act.store,
                             avg.zeta  = if(MH.step)             ifelse(zeta.tune, mean(avgzeta), zeta),
                             time      = init.time)
      return(returns)
  }
