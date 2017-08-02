############################################################################
### Gibbs Sampler for Infinite DP Mixtures of Infinite Factor Analysers ####
############################################################################

# Gibbs Sampler Function
  .gibbs_IMIFA       <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, learn.alpha, mu, sw, uni.type,
                                 uni.prior, sigma.mu, burnin, thinning, a.hyper, psi.alpha, psi.beta, verbose, trunc.G,
                                 adapt, ind.slice, alpha.d1, discount, alpha.d2, cluster, b0, b1, IM.lab.sw, zeta, tune.zeta,
                                 nu, prop, d.hyper, beta.d1, beta.d2, adaptat, epsilon, learn.d, nuplus1, kappa, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    matrix           <- base::matrix
    total            <- max(iters)
    if(verbose)         pb    <- txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
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
    ll.store         <- rep(0L, n.store)
    Q.star           <- Q
    Qs               <- rep(Q,  trunc.G)
    Q.store          <- matrix(0L, nrow=trunc.G, ncol=n.store)
    Q.large          <- Q.big <- Q.bigs <- FALSE
    acc1             <- acc2  <- FALSE
    err.z            <- z.err <- FALSE
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
    MH.step          <- any(discount  > 0, learn.d) && learn.alpha
    if(MH.step)     {
      a.rates        <- rep(0L, total)
    } else a.rates   <- 1
    if(IM.lab.sw)   {
      lab.rate       <- matrix(0L, nrow=2, ncol=total)
    }
    d.count          <- 0
    avgzeta          <- rep(zeta, 100)
    heat             <- tune.zeta$heat
    lambda           <- tune.zeta$lambda
    target           <- tune.zeta$target
    zeta.tune        <- heat > 0 && tune.zeta$do
    mu.sigma         <- 1/sigma.mu
    sig.mu.sqrt      <- sqrt(sigma.mu)
    z                <- cluster$z
    pi.alpha         <- cluster$pi.alpha
    one.uni          <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv     <- switch(uni.type,   unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                           constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip      <- switch(uni.prior,  unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni)) {
      uni.shape      <- switch(uni.type,   constrained=N/2 + psi.alpha, single=(N * P)/2 + psi.alpha)
      V              <- switch(uni.type,   constrained=P, single=1)
    }
    psi.beta         <- switch(uni.prior,  isotropic=as.vector(unique(round(psi.beta, min(nchar(psi.beta))))), psi.beta)
    pi.prop          <- cluster$pi.prop
    nn               <- tabulate(z, nbins=trunc.G)
    mu.tmp           <- vapply(seq_len(trunc.G - G), function(g) .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))
    mu               <- cbind(mu, if(uni) t(mu.tmp) else mu.tmp)
    eta              <- .sim_eta_p(N=N, Q=Q)
    phi              <- lapply(Ts, function(t) .sim_phi_p(Q=Q, P=P, nu=nu, plus1=nuplus1))
    delta            <- lapply(Ts, function(t) c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2)))
    tau              <- lapply(delta, cumprod)
    lmat             <- lapply(Ts, function(t) matrix(vapply(Ps, function(j) .sim_load_ps(Q=Q, phi=phi[[t]][j,], tau=tau[[t]]), numeric(Q)), nrow=P, byrow=TRUE))
    psi.inv          <- vapply(Ts, function(t) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    psi.inv          <- if(uni) t(psi.inv) else psi.inv
    if(Q < min(N - 1, Ledermann(P)))     {
      for(g in which(nn > P)) {
        fact         <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]        <- fact$scores
          lmat[[g]]           <- unclass(fact$loadings)
          psi.inv[,g]         <- 1/fact$uniquenesses
        }
      }
    } else {
      psi.tmp        <- psi.inv
      if(isTRUE(one.uni)) {
        psi.inv[,Gs] <- 1/switch(uni.type, constrained=Rfast::colVars(data), exp(mean(log(Rfast::colVars(data)))))
      } else   {
        psi.inv[,Gs] <- 1/vapply(Gs, function(g) if(nn[g] > 1) switch(uni.type, unconstrained=Rfast::colVars(data[z == g,, drop=FALSE]), rep(exp(mean(log(Rfast::colVars(data[z == g,, drop=FALSE])))), P)) else psi.tmp[,g], numeric(P))
      }
      inf.ind        <- is.infinite(psi.inv)
      psi.inv[inf.ind]        <- psi.tmp[inf.ind]
    }
    index            <- order(pi.prop, decreasing=TRUE)
    pi.prop          <- pi.prop[index]
    mu               <- mu[,index, drop=FALSE]
    phi              <- phi[index]
    delta            <- delta[index]
    tau              <- tau[index]
    lmat             <- lmat[index]
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
    slice.logs       <- c(- Inf, 0L)
    if(burnin         < 1)  {
      if(sw["mu.sw"])   mu.store[,,1]    <- mu
      if(sw["s.sw"])    eta.store[,,1]   <- eta
      if(sw["l.sw"])    load.store[,,,1] <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
      if(sw["psi.sw"])  psi.store[,,1]   <- 1/psi.inv
      if(sw["pi.sw"])   pi.store[,1]     <- pi.prop
      z.store[1,]             <- z
      Q0                      <- Qs > 0
      sigma                   <- if(uni) lapply(Gs, function(g) as.matrix(1/psi.inv[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0)) else lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g]))
      log.probs               <- if(uni) vapply(Gs, function(g) dnorm(data, mu[,g], sqrt(sigma[[g]]), log=TRUE) + log(pi.prop[g]), numeric(N)) else vapply(Gs, function(g, Q=Q0[g]) { sigma <- if(Q) sigma[[g]] else sqrt(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q) + log(pi.prop[g]) }, numeric(N))
      ll.store[1]             <- sum(rowLogSumExps(log.probs))
      Q.store[,1]             <- Qs
      G.store[1]              <- G.non
      act.store[1]            <- G
      if(learn.alpha) {
        alpha.store[1]        <- pi.alpha
      }
      if(learn.d)     {
        d.store[1]            <- discount
      }
    }
    init.time        <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose     && iter   < burnin) setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

    # Mixing Proportions
      Vs             <- .sim_vs_inf(alpha=pi.alpha, nn=nn[Gs], N=N, discount=discount, len=G, lseq=Gs)
      pi.prop        <- .sim_pi_inf(Vs, len=G)

    # Re-ordering & Slice Sampler
      index          <- order(pi.prop[Gs], decreasing=TRUE)
      prev.prod      <- pi.prop[G]  * (1/Vs[G] - 1)
      pi.prop[Gs]    <- pi.prop[index]
      Vs[Gs]         <- Vs[index]
      mu[,Gs]        <- mu[,index, drop=FALSE]
      phi[Gs]        <- phi[index]
      delta[Gs]      <- delta[index]
      tau[Gs]        <- tau[index]
      lmat[Gs]       <- lmat[index]
      Qs[Gs]         <- Qs[index]
      psi.inv[,Gs]   <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.integer(levels(z))[z]
      if(!ind.slice) {
        ksi          <- pi.prop
        log.ksi      <- log(ksi)
      }
      u.slice        <- runif(N, 0, ksi[z])
      min.u          <- min(u.slice)
      G.old          <- G
      if(ind.slice)   {
        G.new        <- sum(min.u   < ksi)
        G.trunc      <- G.new < G.old
        while(G  < G.new && (Vs[G] != 1 || pi.prop[G] != 0)) {
          newVs      <- .sim_vs_inf(alpha=pi.alpha, discount=discount, len=1, lseq=G + 1)
          Vs         <- c(Vs,      newVs)
          pi.prop    <- c(pi.prop, newVs * prev.prod)
          G          <- G + 1
          prev.prod  <- pi.prop[G]  * (1/Vs[G] - 1)
        }
        G            <- ifelse(!G.trunc, G, G.new)
        Gs           <- seq_len(G)
        if(G.trunc)   {
          pi.prop    <- pi.prop[Gs]
          Vs         <- Vs[Gs]
        }
      } else     {
        cum.pi       <- sum(pi.prop)
        u.max        <- 1 - min.u
        G.trunc      <- cum.pi > u.max
        while(cum.pi  < u.max && trunc.G > G && (pi.prop[G] != 0 || Vs[G] != 1)) {
          newVs      <- .sim_vs_inf(alpha=pi.alpha, discount=discount, len=1, lseq=G + 1)
          Vs         <- c(Vs,      newVs)
          newPis     <- newVs  *   prev.prod
          pi.prop    <- c(pi.prop, newPis)
          cum.pi     <- cum.pi +   newPis
          ksi        <- pi.prop
          log.ksi    <- c(log.ksi, log(newPis))
          G          <- G + 1
          prev.prod  <- pi.prop[G] * (1/Vs[G] - 1)
        }
        G            <- ifelse(G.trunc, which.max(cumsum(pi.prop) > u.max), G)
        Gs           <- seq_len(G)
        if(G.trunc) {
          pi.prop    <- ksi   <-   pi.prop[Gs]
          log.ksi    <- log.ksi[Gs]
          Vs         <- Vs[Gs]
        }
      }

    # Cluster Labels
      Qgs            <- Qs[Gs]
      Q0             <- Qgs  > 0
      Q1             <- Qgs == 1
      if(G > 1)     {
        psi          <- 1/psi.inv
        sigma        <- if(uni) lapply(Gs, function(g) as.matrix(psi[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0)) else lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
        if(uni) {
          log.probs  <- vapply(Gs, function(g) dnorm(data, mu[,g], sqrt(sigma[[g]]), log=TRUE), numeric(N))
        } else  {
          log.check  <- capture.output(log.probs <- try(vapply(Gs, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q), numeric(N)), silent=TRUE))
        }
        if(inherits(log.probs, "try-error")) {
          log.probs  <- vapply(Gs, function(g, Q=Q0[g]) { sigma <- if(Q) sigma[[g]] else sqrt(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q) }, numeric(N))
        }
        if(ind.slice) {
          log.pixi   <- log(pi.prop) - log(ksi[Gs])
          log.probs  <- vapply(Gs, function(g) slice.logs[1 + (u.slice < ksi[g])] + log.pixi[g], numeric(N)) + log.probs
        } else   {
          log.probs  <- vapply(Gs, function(g) slice.logs[1 + (u.slice < ksi[g])], numeric(N)) + log.probs
        }
        z            <- gumbel_max(probs=log.probs, slice=TRUE)
      } else   {
        z            <- rep(1L, N)
      }
      nn             <- tabulate(z, nbins=trunc.G)
      nn0            <- nn > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)
      dat.g          <- lapply(Gs, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data         <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta          <- .empty_mat(N)
        eta.tmp      <- lapply(Gs, function(g) eta[z == g,, drop=FALSE])
        lmat[Gs]     <- replicate(G, .empty_mat(P))
      } else {
        eta.tmp      <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else matrix(0, nrow=ifelse(Q0[g], 0, nn[g]), ncol=Qs[g]))
        EtE          <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[Gs]     <- lapply(Gs, function(g) matrix(if(all(nn0[g], Q0[g])) vapply(Ps, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j],
                               Q1=Q1[g], EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]]), numeric(Qs[g])) else
                               vapply(Ps, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]]), numeric(Qs[g])), nrow=P, byrow=TRUE))
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat        <- lapply(Gs, function(g) { S <- c.data[[g]] - if(Q0[g]) tcrossprod(eta.tmp[[g]], lmat[[g]]) else 0; S * S } )
        psi.inv[,]   <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,Gs] <- vapply(Gs, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                        P=P, Q0=Q0[g], eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      }

    # Means
      sum.data       <- vapply(dat.g, colSums, numeric(P))
      sum.data       <- if(uni) t(sum.data) else sum.data
      sum.eta        <- lapply(eta.tmp, colSums)
      mu[,Gs]        <- vapply(Gs, function(g) if(nn0[g]) .sim_mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                               sum.data=sum.data[,g], lmat=lmat[[g]], mu.zero=mu.zero) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Shrinkage
      if(all(Q0))     {
        load.2       <- lapply(lmat[Gs], .power2)
        phi[Gs]      <- lapply(Gs, function(g) if(nn0[g]) .sim_phi(Q=Qs[g], P=P, nu=nu, plus1=nuplus1, tau=tau[[g]],
                               load.2=load.2[[g]]) else .sim_phi_p(Q=Qs[g], P=P, nu=nu, plus1=nuplus1))

        sum.terms    <- lapply(Gs, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
        for(g in Gs)  {
          Qg         <- Qs[g]
          Q1g        <- Q1[g]
          if(nn0[g])  {
            for(k in seq_len(Qg)) {
              delta[[g]][k]   <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P=P,
                                 Q=Qg, k=k, sum.term.kq=sum.terms[[g]][k:Qg]) else .sim_delta1(Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]],
                                 alpha.d1=ifelse(Q1g, alpha.d2, alpha.d1), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1])
              tau[[g]]        <- cumprod(delta[[g]])
          }
          } else {
            for(k in seq_len(Qg)) {
              delta[[g]][k]   <- if(k > 1) .sim_delta_p(alpha=alpha.d2, beta=beta.d2) else .sim_delta_p(alpha=ifelse(Q1g, alpha.d2, alpha.d1), beta=ifelse(Q1g, beta.d2, beta.d1))
              tau[[g]]        <- cumprod(delta[[g]])
            }
          }
        }
      }

    # Adaptation
      if(all(adapt, iter > adaptat)) {
        if(runif(1)   < ifelse(iter  < burnin, 0.5, exp(-b0 - b1 * (iter - adaptat)))) {
          colvec     <- lapply(nn.ind, function(g) (if(Q0[g]) colSums(abs(lmat[[g]])   < epsilon)/P else 0) >= prop)
          nonred     <- lapply(colvec, .which0)
          numred     <- lengths(colvec) - lengths(nonred)
          notred     <- numred == 0
          ng.ind     <- seq_along(nn.ind)
          Qs.old     <- Qs[nn0]
          Qs[nn0]    <- pmax(0, vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1 else Qs.old[h] - numred[h], numeric(1L)))
          Q.big      <- Qs[nn0] > Q.star
          if((Q.bigs <- any(Q.big))) {
            notred   <- notred & !Q.big
            Qs[nn0][Q.big]    <- Q.star
          }
          phi[nn0]   <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])],  rgamma(n=P, shape=nu + nuplus1, rate=nu)) else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0] <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])],     rgamma(n=1, shape=alpha.d2, rate=beta.d2)) else delta[[g]][nonred[[h]]])
          tau[nn0]   <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]  <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp       <- Qs[!nn0]
          Qpop       <- Qs[nn0]
          Qmax       <- ifelse(all(Q.big), max(Qpop), max(Qpop[!Q.big]))
          Qmaxold    <- max(Qs.old)
          store.eta  <- all(sw["s.sw"], storage)
          if(any(!nn0)    && Qmax  !=  max(Qemp)) {
            Qmaxseq  <- seq_len(Qmax)
            for(t in Ts[!nn0][Qemp !=  Qmax])     {
              Qt     <- Qs[t]
              if(Qt   > Qmax)  {
                phi[[t]]      <- phi[[t]][,Qmaxseq,  drop=FALSE]
                delta[[t]]    <- delta[[t]][Qmaxseq]
                tau[[t]]      <- tau[[t]][Qmaxseq]
                lmat[[t]]     <- lmat[[t]][,Qmaxseq, drop=FALSE]
              } else  {
                while(Qt  != Qmax)  {
                 phi[[t]]     <- cbind(phi[[t]],     rgamma(n=P, shape=nu + nuplus1, rate=nu))
                 delta[[t]]   <- c(delta[[t]],       rgamma(n=1, shape=alpha.d2, rate=beta.d2))
                 tau[[t]]     <- cumprod(delta[[t]])
                 if(store.eta && t %in% Gs)   {
                 eta.tmp[[t]] <- cbind(eta.tmp[[t]], matrix(0, nr=0, nc=1))
                 }
                 Qt  <- Qt + 1
                 lmat[[t]]    <- cbind(lmat[[t]],    rnorm(n=P, mean=0, sd=sqrt(1/(phi[[t]][,Qt] * tau[[t]][Qt]))))
                }
              }
            }
            Qs[Qmax  != Qs & !nn0] <-  Qmax
          }
          if(store.eta)    {
            eta.tmp  <- lapply(Gs,     function(g) if(nn0[g] && Qs[g] > Qs.old[which(nn.ind == g)]) cbind(eta.tmp[[g]], rnorm(nn[g])) else eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE])
          }
        }
      }

    # Alpha
      if(learn.alpha)      {
        non.zero.d   <- discount != 0
        if(isTRUE(non.zero.d))  {
          MH.alpha   <- .sim_alpha_m(alpha=pi.alpha, discount=discount, alpha.shape=alpha.shape, alpha.rate=alpha.rate, N=N, G=G.non, zeta=zeta)
          pi.alpha   <- MH.alpha$alpha
          a.rate     <- MH.alpha$rate
          if(isTRUE(zeta.tune)) {
            d.count  <- d.count + non.zero.d
            if(iter   > 100)    {
              zeta   <- .tune_zeta(zeta=zeta, time=d.count, l.rate=MH.alpha$l.prob, heat=heat, target=target, lambda=lambda)
            }
            avgzeta  <- c(avgzeta, zeta)
          }
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
          if((acc1   <- move1$rate1)) {
            sw1      <- move1$sw
            sw1x     <- c(sw1[2], sw1[1])
            nn[sw1]  <- nn[sw1x]
            nn0[sw1] <- nn0[sw1x]
            nn.ind   <- which(nn0)
            Vs[sw1]  <- Vs[sw1x]
            mu[,sw1] <- mu[,sw1x, drop=FALSE]
            phi[sw1] <- phi[sw1x]
            tau[sw1] <- tau[sw1x]
            Qs[sw1]  <- Qs[sw1x]
            lmat[sw1]         <- lmat[sw1x]
            delta[sw1]        <- delta[sw1x]
            psi.inv[,sw1]     <- psi.inv[,sw1x, drop=FALSE]
            pi.prop[sw1]      <- pi.prop[sw1x]
            zsw1     <- z == sw1[1]
            z[z == sw1[2]]    <- sw1[1]
            z[zsw1]  <- sw1[2]
            if(all(sw["s.sw"], storage)) {
              eta.tmp[sw1]    <- eta.tmp[sw1x]
              dat.g[sw1]      <- dat.g[sw1x]
              Q0[sw1]         <- Q0[sw1x]
            }
          }
        } else  acc1 <- FALSE
        if(G     > 1) {
          move2      <- .lab_move2(G=G, Vs=Vs, nn=nn)
          if((acc2   <- move2$rate2)) {
            sw2      <- move2$sw
            sw2x     <- c(sw2[2], sw2[1])
            nn[sw2]  <- nn[sw2x]
            nn0[sw2] <- nn0[sw2x]
            nn.ind   <- which(nn0)
            mu[,sw2] <- mu[,sw2x, drop=FALSE]
            phi[sw2] <- phi[sw2x]
            tau[sw2] <- tau[sw2x]
            Qs[sw2]  <- Qs[sw2x]
            lmat[sw2]         <- lmat[sw2x]
            delta[sw2]        <- delta[sw2x]
            psi.inv[,sw2]     <- psi.inv[,sw2x, drop=FALSE]
            pi.prop[sw2]      <- pi.prop[sw2x]
            zsw2     <- z == sw2[1]
            z[z == sw2[2]]    <- sw2[1]
            z[zsw2]  <- sw2[2]
            if(all(sw["s.sw"], storage)) {
              eta.tmp[sw2]    <- eta.tmp[sw2x]
              dat.g[sw2]      <- dat.g[sw2x]
              Q0[sw2]         <- Q0[sw2x]
            }
          }
        } else  acc2 <- FALSE
      }

      if(Q.bigs && !Q.large   && iter > burnin) {         warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large      <- TRUE
      }
      if(z.err  && !err.z) {                              warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z        <- TRUE
      }
      if(MH.step)        a.rates[iter]                 <- a.rate
      if(learn.d)        d.rates[iter]                 <- d.rate
      if(IM.lab.sw)      lab.rate[,iter]               <- c(acc1, acc2)
      if(storage)    {
        if(verbose)      setTxtProgressBar(pb, iter)
        new.it       <-  which(iters == iter)
        if(sw["mu.sw"])  mu.store[,,new.it]            <- mu
        if(all(sw["s.sw"],
           any(Q0)))    {
          max.Q      <-  max(Qs)
          eta.tmp    <-  if(length(unique(Qs)) != 1)      lapply(Gs,       function(g) cbind(eta.tmp[[g]], matrix(0, nrow=nn[g], ncol=max.Q - Qs[g]))) else eta.tmp
          q0ng       <-  (!Q0  | Qs[Gs] == 0)   & nn0[Gs]
          if(any(q0ng)) {
            eta.tmp[q0ng]     <-                          lapply(Gs[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                   eta.store[,seq_len(max.Q),new.it]   <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])  {
          for(g in Gs)  {
            Qseqg    <- seq_len(Qs[g])
            if(Q0[g])   load.store[,Qseqg,g,new.it]    <- lmat[[g]]
          }
        }
        if(sw["psi.sw"]) psi.store[,,new.it]           <- 1/psi.inv
        if(sw["pi.sw"])  pi.store[Gs,new.it]           <- pi.prop
        if(learn.alpha)  alpha.store[new.it]           <- pi.alpha
        if(learn.d)      d.store[new.it]               <- discount
                         z.store[new.it,]              <- as.integer(z)
                         ll.store[new.it]              <- if(G > 1) sum(rowLogSumExps(log.probs)) else sum(dmvn(X=data, mu=mu[,nn.ind], sigma=tcrossprod(lmat[[nn.ind]]) + diag(psi.store[,nn.ind,new.it]), log=TRUE))
                         Q.store[,new.it]              <- as.integer(Qs)
                         G.store[new.it]               <- as.integer(G.non)
                         act.store[new.it]             <- as.integer(G)
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
    returns          <- list(mu        = if(sw["mu.sw"])   tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                             eta       = if(sw["s.sw"])    tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                             load      = if(sw["l.sw"])    tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                             psi       = if(sw["psi.sw"])  tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                             pi.prop   = if(sw["pi.sw"])   pi.store,
                             alpha     = if(learn.alpha)   alpha.store,
                             discount  = if(learn.d) {     if(sum(d.store == 0)/n.store > 0.5) as.simple_triplet_matrix(d.store) else d.store },
                             a.rate    = ifelse(MH.step,   mean(a.rates), a.rates),
                             d.rate    = ifelse(learn.d,   mean(d.rates), d.rates),
                             lab.rate  = if(IM.lab.sw)     setNames(rowmeans(lab.rate), c("Move1", "Move2")),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             Q.store   = tryCatch(Q.store[Gmax,, drop=FALSE],          error=function(e) Q.store),
                             G.store   = G.store,
                             act.store = act.store,
                             avg.zeta  = if(MH.step)       ifelse(zeta.tune, mean(avgzeta), zeta),
                             time      = init.time)
    attr(returns, "Q.big")    <- Q.large
    return(returns)
  }
