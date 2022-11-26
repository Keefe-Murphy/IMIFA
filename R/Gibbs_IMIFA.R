############################################################################
### Gibbs Sampler for Infinite DP Mixtures of Infinite Factor Analysers ####
############################################################################

# Gibbs Sampler Function
  .gibbs_IMIFA       <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, learn.alpha, mu, sw, uni.type, uni.prior,
                                 sigma.mu, burnin, thinning, a.hyper, psi.alpha, psi.beta, verbose, trunc.G, adapt, ind.slice, discount,
                                 alpha.d1, alpha.d2, cluster, b0, b1, IM.lab.sw, zeta, tune.zeta, rho1, rho2, nu1, nu2, truncated, cluster.shrink,
                                 prop, d.hyper, beta.d1, beta.d2, start.AGS, stop.AGS, epsilon, learn.d, kappa, forceQg, thresh, exchange, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    sq_mat           <- if(P   > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix           <- base::matrix
    total            <- max(iters)
    if(verbose)         pb    <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    AGS.burn         <- total/5L
    Gs               <- seq_len(G)
    Ts               <- seq_len(trunc.G)
    Ps               <- seq_len(P)
    Ns               <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
    uni              <- P == 1
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
    Q.store          <- matrix(0L, nrow=trunc.G, ncol=n.store)
    Q.large          <- Q.big <- Q.bigs <- FALSE
    acc1             <- acc2  <- FALSE
    err.z            <- z.err <- FALSE
    act.store        <- G.store
    pi.alpha         <- cluster$pi.alpha
    nu1.5            <- nu1 + 0.5
    P.5              <- P/2
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
    z                <- cluster$z
    nn               <- tabulate(z, nbins=trunc.G)
    nn0              <- nn  > 0
    nn.ind           <- which(nn > 0)
    G.non            <- length(nn.ind)
    Q.star           <- Q
    Qs               <- rep(Q, trunc.G)
    Qs               <- if(forceQg) pmin(Qs, replace(nn, !nn0, Inf) - 1L) else Qs
    Q0               <- Qs  > 0
    Q1               <- Qs == 1
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
    Qmax             <- ifelse(forceQg, max(Qs), Q)
    Qmaxseq          <- seq_len(Qmax)
    eta              <- .sim_eta_p(N=N, Q=Qmax)
    phi              <- if(forceQg) lapply(Ts, function(t) .sim_phi_p(Q=Qs[t], P=P, nu1=nu1, nu2=nu2)) else replicate(trunc.G, .sim_phi_p(Q=Q, P=P, nu1=nu1, nu2=nu2), simplify=FALSE)
    if(isTRUE(truncated))   {
      .sim_deltak    <- .sim_deltaKT
      .rdelta        <- rltrgamma
      delta          <- if(forceQg) lapply(Ts, function(t) c(if(Qs[t] > 0) .sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_deltaPT(Q=Qs[t], alpha=alpha.d2, beta=beta.d2))) else replicate(trunc.G, list(c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_deltaPT(Q=Q, alpha=alpha.d2, beta=beta.d2))))
      .sim_delta_p   <- .sim_deltaPT
    } else              {
      .rdelta        <- stats::rgamma
      delta          <- if(forceQg) lapply(Ts, function(t) c(if(Qs[t] > 0) .sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Qs[t], alpha=alpha.d2, beta=beta.d2))) else replicate(trunc.G, list(c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2))))
    }
    tau              <- lapply(delta, cumprod)
    if(cluster.shrink)  {
      sig.store      <- matrix(0L, nrow=trunc.G, ncol=n.store)
      MGPsig         <- .sim_sigma_p(G=trunc.G, rho1=rho1, rho2=rho2)
    } else MGPsig    <- rep(1L, trunc.G)
    lmat             <- lapply(Ts, function(t) matrix(vapply(Ps, function(j) .sim_load_ps(Q=Qs[t], phi=phi[[t]][j,], tau=tau[[t]], sigma=MGPsig[t]), numeric(Qs[t])), nrow=P, byrow=TRUE))
    if(isTRUE(one.uni)) {
      psi.inv        <- matrix(.sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), nrow=P, ncol=trunc.G)
    } else psi.inv   <- replicate(trunc.G, .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), simplify="array")
    psi.inv          <- if(uni) t(psi.inv) else psi.inv
    if(isTRUE(one.uni)) {
      psi.inv[]      <- 1/switch(EXPR=uni.type, constrained=colVars(data), max(colVars(data)))
    } else   {
      tmp.psi        <- (nn[nn0] - 1L)/pmax(rowsum(data^2, z) - rowsum(data, z)^2/nn[nn0], 0L)
      tmp.psi        <- switch(EXPR=uni.type, unconstrained=t(tmp.psi), matrix(Rfast::rowMaxs(tmp.psi, value=TRUE), nrow=P, ncol=G, byrow=TRUE))
      psi.inv[,nn     > 1]    <- tmp.psi[!is.nan(tmp.psi)]
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
      if(verbose     && iter     < burnin) utils::setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

      # Adaptation
      if(adapt       && all(iter >= start.AGS, iter < stop.AGS))  {
        if(stats::runif(1) < ifelse(iter < AGS.burn, 0.5, exp(-b0 - b1 * (iter - start.AGS))))    {
          colvec     <- lapply(nn.ind, function(g) if(Q0[g]) (colSums2(abs(lmat[[g]])   < epsilon)/P) >= prop else stats::runif(1) <= prop)
          nonred     <- lapply(colvec, .which0)
          numred     <- lengths(colvec)  - lengths(nonred)
          notred     <- numred == 0
          ng.ind     <- seq_along(nn.ind)
          Qs.old     <- Qs[nn0]
          Qs[nn0]    <- pmax.int(0L, vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1L         else Qs.old[h] - numred[h], numeric(1L)))
          star.Q     <- if(forceQg) pmin(Q.star, nn[nn0] - 1L) else Q.star
          Q.big      <- Qs[nn0] > star.Q
          if((Q.bigs <- any(Q.big)))  {
            notred   <- notred & !Q.big
            Qs[nn0][Q.big]    <- if(forceQg) star.Q[Q.big]     else Q.star
            if(forceQg)        {
              for(qb in which(Q.big)) {
                nonred[[qb]]  <- seq_len(star.Q[qb])
              }
            }
          }
          phi[nn0]   <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])],  .rgamma0(n=P, shape=nu1,      rate=nu2))     else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0] <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])],     .rdelta(n=1,  shape=alpha.d2, rate=beta.d2)) else delta[[g]][nonred[[h]]])
          tau[nn0]   <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]  <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]] * MGPsig[g]))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp       <- Qs[!nn0]
          Qmax       <- max(Qs[nn0])
          Qmaxseq    <- seq_len(Qmax)
          if(any(!nn0)     && max(Qemp) != Qmax)  {
            for(t in Ts[!nn0][Qemp      != Qmax]) {
              Qt     <- Qs[t]
              if(Qt   > Qmax)  {
                phi[[t]]      <- phi[[t]][,Qmaxseq,  drop=FALSE]
                delta[[t]]    <- delta[[t]][Qmaxseq]
                tau[[t]]      <- tau[[t]][Qmaxseq]
                lmat[[t]]     <- lmat[[t]][,Qmaxseq, drop=FALSE]
              } else  {
                while(Qt   != Qmax)  {
                 phi[[t]]     <- cbind(phi[[t]],     .rgamma0(n=P, shape=nu1,       rate=nu2))
                 delta[[t]]   <- c(delta[[t]],       .rdelta(n=1,  shape=alpha.d2,  rate=beta.d2))
                 tau[[t]]     <- cumprod(delta[[t]])
                 Qt  <- Qt  + 1L
                 lmat[[t]]    <- cbind(lmat[[t]],    stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[[t]][,Qt] * tau[[t]][Qt] * MGPsig[t])))
                }
              }
            }
            Qs[Qmax  != Qs  & !nn0] <- Qmax
          }
          Q0         <- Qs  > 0
          Q1         <- Qs == 1
        }
      }

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
      phi[Gs]        <- phi[index]
      delta[Gs]      <- delta[index]
      tau[Gs]        <- tau[index]
      MGPsig[Gs]     <- MGPsig[index]
      lmat[Gs]       <- lmat[index]
      Qs[Gs]         <- Qs[index]
      Q0[Gs]         <- Q0[index]
      Q1[Gs]         <- Q1[index]
      psi.inv[,Gs]   <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.integer(levels(z))[z]
      nn[Gs]         <- nn[index]
      nn0[Gs]        <- nn0[index]

    # Scores & Loadings
      dat.g          <- lapply(Gs, function(g) data[z == g,, drop=FALSE])
      c.data         <- lapply(Gs, function(g) sweep(dat.g[[g]], 2L, mu[,g], FUN="-", check.margin=FALSE))
      n.eta          <- nn
      n0q0           <- nn0  & Q0
      q0ng           <- (!Q0 | Q1) & n.eta > 0
      if(!any(Q0))    {
        eta          <- .empty_mat(nr=N)
        eta.tmp      <- lapply(Gs, function(g) eta[z == g,,  drop=FALSE])
        lmat[Gs]     <- replicate(G, .empty_mat(nr=P))
      } else {
        eta.tmp      <- lapply(Gs, function(g) if(n0q0[g]) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0L, nrow=ifelse(Q0[g], 0L, nn[g]), ncol=Qs[g]))
        EtE          <- lapply(Gs, function(g) if(n0q0[g]) crossprod(eta.tmp[[g]]))
        lmat[Gs]     <- lapply(Gs, function(g) matrix(if(n0q0[g]) vapply(Ps, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j], Q1=Q1[g],
                               EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g]))   else
                               vapply(Ps, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g])), nrow=P, byrow=TRUE))
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat        <- lapply(Gs, function(g) { S <- c.data[[g]] - if(Q0[g]) tcrossprod(eta.tmp[[g]], lmat[[g]]) else 0L; S^2 } )
        psi.inv[,]   <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,Gs] <- vapply(Gs, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                        P=P, Q0=Q0[g], eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      }

    # Means
      sum.data       <- vapply(dat.g, colSums2, numeric(P))
      sum.data       <- if(uni) t(sum.data) else sum.data
      sum.eta        <- lapply(eta.tmp, colSums2)
      mu[,Gs]        <- vapply(Gs, function(g) if(nn0[g]) .sim_mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                               sum.data=sum.data[,g], lmat=lmat[[g]], mu.prior=mu.prior) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Shrinkage
      if(any(Q0))     {
        load.2       <- lapply(lmat[Gs], "^", 2)
        phi[Gs]      <- lapply(Gs, function(g) if(n0q0[g]) .sim_phi(Q=Qs[g], P=P, nu1.5=nu1.5, nu2=nu2, tau=tau[[g]],
                        load.2=load.2[[g]], sigma=MGPsig[g]) else .sim_phi_p(Q=Qs[g], P=P, nu1=nu1, nu2=nu2))
        sum.terms    <- lapply(Gs, function(g) if(n0q0[g]) colSums2(phi[[g]] * load.2[[g]]))
        for(g in Gs)  {
          Qg         <- Qs[g]
          Q1g        <- Q1[g]
          if(n0q0[g]) {
            for(k in seq_len(Qg)) {
              delta[[g]][k]   <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P.5=P.5, Q=Qg,
                                 k=k, sum.term.kq=sum.terms[[g]][k:Qg], sigma=MGPsig[g]) else .sim_delta1(Q=Qg, P.5=P.5, tau=tau[[g]], sum.term=sum.terms[[g]],
                                 alpha.d1=ifelse(Q1g, alpha.d2, alpha.d1), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1L], sigma=MGPsig[g])
              tau[[g]]        <- cumprod(delta[[g]])
            }
          } else {
            if(Q0[g])  {
              delta[[g]]      <-  c(stats::rgamma(n=1, shape=ifelse(Q1g, alpha.d2, alpha.d1), rate=ifelse(Q1g, beta.d2, beta.d1)),
                                    .sim_delta_p(Q=Qg, alpha=alpha.d2, beta=beta.d2))
              tau[[g]]        <- cumprod(delta[[g]])
            }
          }
        }
        if(cluster.shrink)     {
          nnX                 <- n0q0[Gs]
          n0Gq                <- which(nnX)
          nGq0                <- length(n0Gq)
          MGPsig[n0Gq]        <- .sim_sigma(G=nGq0, P.5=P.5, Qs=Qs[n0Gq], rho1=rho1, rho2=rho2, sum.terms=sum.terms[n0Gq], tau=tau[n0Gq])
          MGPsig[which(!nnX)] <- .sim_sigma_p(G=G - nGq0, rho1=rho1, rho2=rho2)
        }
      }

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
        if(G.trunc) {
          pi.prop    <- ksi   <-   pi.prop[Gs]
          Vs         <- Vs[Gs]
        }
      }
      if(thresh)      {
        ksi          <- pmax(pi.prop, TRX)
        log.ksi      <- log(ksi)
      }

    # Cluster Labels
      if(G > 1)       {
        psi          <- 1/psi.inv
        sigma        <- if(uni) lapply(Gs, function(g) as.matrix(psi[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0L)) else lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
        if(ind.slice) {
          log.pixi   <- log(pi.prop) - log.ksi[Gs]
          log.probs  <- vapply(Gs, function(g) { slinf[u.slice < ksi[g]] <- log.pixi[g]; slinf }, numeric(N))
        } else  {
          log.probs  <- vapply(Gs, function(g) slinf[1L + (u.slice < ksi[g])], numeric(N))
        }
        fin.probs    <- is.finite(log.probs)
        if(uni) {
          log.probs  <- vapply(Gs, function(g, LP=log.probs[,g], FP=fin.probs[,g]) { LP[FP] <- stats::dnorm(data[FP,], mu[,g], sq_mat(sigma[[g]]), log=TRUE) + LP[FP]; LP }, numeric(N))
        } else  {
          probs.log  <- log.probs
          log.probs  <- try(vapply(Gs, function(g, LP=log.probs[,g], FP=fin.probs[,g]) { LP[FP] <- dmvn(data[FP,], mu[,g], sigma[[g]], log=TRUE,  isChol=!Q) + LP[FP]; LP }, numeric(N)), silent=TRUE)
          if(z.err   <- inherits(log.probs, "try-error")) {
           log.probs <- vapply(Gs, function(g, LP=probs.log[,g], FP=fin.probs[,g], Q=Q0[g]) { sigma <- if(Q) is.posi_def(sigma[[g]], make=TRUE)$X.new else sq_mat(sigma[[g]]); LP[FP] <- dmvn(data[FP,], mu[,g], sigma, log=TRUE, isChol=!Q) + LP[FP]; LP }, numeric(N))
          }
          rm(probs.log)
        }
        z            <- gumbel_max(probs=log.probs, slice=TRUE)
      } else    {
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
            phi[sw1] <- phi[sw1x]
            tau[sw1] <- tau[sw1x]
            Qs[sw1]  <- Qs[sw1x]
            MGPsig[sw1]       <- MGPsig[sw1x]
            lmat[sw1]         <- lmat[sw1x]
            delta[sw1]        <- delta[sw1x]
            psi.inv[,sw1]     <- psi.inv[,sw1x, drop=FALSE]
            pi.prop[sw1]      <- pi.prop[sw1x]
            zsw1     <- z == sw1[1L]
            z[z == sw1[2L]]   <- sw1[1L]
            z[zsw1]  <- sw1[2L]
            Q0[sw1]  <- Q0[sw1x]
            Q1[sw1]  <- Q1[sw1x]
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
            phi[sw2] <- phi[sw2x]
            tau[sw2] <- tau[sw2x]
            Qs[sw2]  <- Qs[sw2x]
            MGPsig[sw2]       <- MGPsig[sw2x]
            lmat[sw2]         <- lmat[sw2x]
            delta[sw2]        <- delta[sw2x]
            psi.inv[,sw2]     <- psi.inv[,sw2x, drop=FALSE]
            pi.prop[sw2]      <- pi.prop[sw2x]
            zsw2     <- z == sw2[1L]
            z[z == sw2[2L]]   <- sw2[1L]
            z[zsw2]  <- sw2[2L]
            Q0[sw2]  <- Q0[sw2x]
            Q1[sw2]  <- Q1[sw2x]
          }
        } else  acc2 <- FALSE
      }

      if(Q.bigs && !Q.large   && iter > burnin) {         cat("\n"); cat("\n"); warning(paste0("\nQ has exceeded initial number of loadings columns", ifelse(forceQg, " (or exceeded the number of observations in one or more clusters)", ""), " since burnin:\nconsider increasing 'range.Q' from ", Q.star, ifelse(forceQg, " or setting 'forceQg' to FALSE\n", "\n")), call.=FALSE)
        Q.large      <- TRUE
      }
      if(z.err  && !err.z) {                              cat("\n"); warning("\nAlgorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices\n", call.=FALSE)
        err.z        <- TRUE
      }
      if(MH.step)        a.rates[iter]                 <- a.rate
      if(learn.d)        d.rates[iter]                 <- d.rate
      if(IM.lab.sw)      lab.rate[,iter]               <- c(acc1, acc2)
      if(storage)    {
        if(verbose)      utils::setTxtProgressBar(pb, iter)
        new.it       <-  which(iters == iter)
        if(sw["mu.sw"])          mu.store[,,new.it]    <-   mu
        if(all(sw["s.sw"],
           any(Q0)))    {
          eta.tmp    <-  if(length(unique(Qs)) != 1)        lapply(seq_len(G.old),       function(g) cbind(eta.tmp[[g]], base::matrix(0L, nrow=n.eta[g], ncol=Qmax - ncol(eta.tmp[[g]])))) else eta.tmp
          if(any(q0ng)) {
            eta.tmp[q0ng]     <-                            lapply(seq_len(G.old)[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                         eta.store[,Qmaxseq,new.it]    <-   do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])  {
          for(g in Gs)  {
            Qseqg    <- seq_len(Qs[g])
            if(Q0[g])   load.store[,Qseqg,g,new.it]    <-   lmat[[g]]
          }
        }
        if(sw["psi.sw"])        psi.store[,,new.it]    <-   1/psi.inv
        if(sw["pi.sw"])         pi.store[Gs,new.it]    <-   pi.prop
        if(learn.alpha)         alpha.store[new.it]    <-   pi.alpha
        if(learn.d)                 d.store[new.it]    <-   discount
        if(cluster.shrink)       sig.store[,new.it]    <-   MGPsig
                                   z.store[new.it,]    <-   as.integer(z)
                                   ll.store[new.it]    <-   if(G > 1) sum(rowLogSumExps(log.probs)) else sum(dmvn(X=data, mu=mu[,nn.ind], sigma=tcrossprod(lmat[[nn.ind]]) + diag(psi.store[,nn.ind,new.it]), log=TRUE))
                                   Q.store[,new.it]    <-   as.integer(Qs)
                                    G.store[new.it]    <-   as.integer(G.non)
                                  act.store[new.it]    <-   as.integer(G)
      }
    }
    if(verbose)         close(pb)
    Gmax             <- seq_len(max(as.integer(z.store)))
    Qmax             <- seq_len(max(Q.store))
    mu.store         <- if(sw["mu.sw"])  tryCatch(mu.store[,Gmax,, drop=FALSE],        error=function(e) mu.store)
    eta.store        <- if(sw["s.sw"])   tryCatch(eta.store[,Qmax,, drop=FALSE],       error=function(e) eta.store)
    load.store       <- if(sw["l.sw"])   tryCatch(load.store[,Qmax,Gmax,, drop=FALSE], error=function(e) load.store)
    psi.store        <- if(sw["psi.sw"]) tryCatch(psi.store[,Gmax,, drop=FALSE],       error=function(e) psi.store)
    pi.store         <- if(sw["pi.sw"])  tryCatch(pi.store[Gmax,, drop=FALSE],         error=function(e) pi.store)
    returns          <- list(mu        = if(sw["mu.sw"])    tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                             eta       = if(sw["s.sw"])     tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                             load      = if(sw["l.sw"])     tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                             psi       = if(sw["psi.sw"])   tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                             pi.prop   = if(sw["pi.sw"])    pi.store,
                             sigma     = if(cluster.shrink) sig.store,
                             alpha     = if(learn.alpha)    alpha.store,
                             discount  = if(learn.d) {      if(sum(d.store == 0)/n.store > 0.5) as.simple_triplet_matrix(d.store) else d.store },
                             a.rate    = ifelse(MH.step,    mean(a.rates), a.rates),
                             d.rate    = ifelse(learn.d,    mean(d.rates), d.rates),
                             lab.rate  = if(IM.lab.sw)      stats::setNames(rowMeans2(lab.rate), c("Move1", "Move2")),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             Q.store   = tryCatch(Q.store[Gmax,, drop=FALSE],          error=function(e) Q.store),
                             G.store   = G.store,
                             act.store = act.store,
                             avg.zeta  = if(MH.step)        ifelse(zeta.tune, mean(avgzeta), zeta),
                             time      = init.time)
    attr(returns, "Q.big")    <- Q.large
      return(returns)
  }
