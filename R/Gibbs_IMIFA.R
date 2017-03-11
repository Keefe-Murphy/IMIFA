############################################################################
### Gibbs Sampler for Infinite DP Mixtures of Infinite Factor Analysers ####
############################################################################

# Gibbs Sampler Function
  .gibbs_IMIFA       <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, learn.alpha, mu, sw, uni.type,
                                 uni.prior, sigma.mu, burnin, thinning, a.hyper, psi.alpha, psi.beta, verbose,
                                 adapt, ind.slice, alpha.d1, discount, alpha.d2, cluster, b0, b1, DP.lab.sw, zeta,
                                 nu, prop, d.hyper, beta.d1, beta.d2, adapt.at, epsilon, learn.d, nuplus1, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    total            <- max(iters)
    if(verbose)         pb    <- txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    trunc.G          <- G
    Gs               <- Ts    <- seq_len(G)
    Ps               <- seq_len(P)
    Ns               <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
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
    z.store          <- matrix(0, nrow=N, ncol=n.store)
    ll.store         <- rep(0, n.store)
    Q.star           <- Q
    Qs               <- rep(Q, trunc.G)
    Q.store          <- matrix(0, nrow=trunc.G, ncol=n.store)
    Q.large          <- Q.big <- Q.bigs <- FALSE
    acc1             <- acc2  <- FALSE
    err.z            <- z.err <- FALSE
    G.store          <- ll.store
    act.store        <- ll.store
    if(learn.alpha) {
      alpha.store    <- ll.store
      alpha.shape    <- a.hyper[1]
      alpha.rate     <- a.hyper[2]
    }
    if(learn.d)     {
      d.store        <- ll.store
      d.shape1       <- d.hyper[1]
      d.shape2       <- d.hyper[2]
      d.rates        <- ll.store
    } else d.rates   <- 1
    MH.step          <- any(discount > 0, learn.d)
    if(MH.step)     {
      a.rates        <- ll.store
    } else a.rates   <- 1
    if(DP.lab.sw) {
      lab.rate       <- matrix(0, nrow=2, ncol=n.store)
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
    eta              <- .sim_eta_p(N=N, Q=Q)
    phi              <- lapply(Ts, function(t) .sim_phi_p(Q=Q, P=P, nu=nu, plus1=nuplus1))
    delta            <- lapply(Ts, function(t) c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2)))
    tau              <- lapply(delta, cumprod)
    lmat             <- lapply(Ts, function(t) matrix(unlist(lapply(Ps, function(j) .sim_load_ps(Q=Q, phi=phi[[t]][j,], tau=tau[[t]])), use.names=FALSE), nrow=P, byrow=TRUE))
    psi.inv          <- vapply(Ts, function(t) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q < .ledermann(N, P))  {
      for(g in which(nn > P)) {
        fact         <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]        <- fact$scores
          lmat[[g]]           <- fact$loadings
          psi.inv[,g]         <- 1/fact$uniquenesses
        }
      }
    } else {
      psi.tmp        <- psi.inv
      psi.inv        <- vapply(Gs, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
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
    z                <- as.numeric(levels(z))[z]
    ksi              <- (1 - rho) * rho^(Ts - 1)
    log.ksi          <- log(ksi)
    slice.logs       <- c(- Inf, 0)
    if(burnin         < 1)  {
      mu.store[,,1]           <- mu
      eta.store[,,1]          <- eta
      load.store[,,,1]        <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
      psi.store[,,1]          <- 1/psi.inv
      pi.store[,1]            <- pi.prop
      z.store[,1]             <- z
      sigma                   <- lapply(Gs, function(g) make.positive.definite(tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])))
      Q0                      <- Qs > 0
      log.probs               <- vapply(Gs, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      ll.store[1]             <- sum(gumbel_max(probs=log.probs, log.like=TRUE)$log.like)
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
      weights        <- .sim_pi_inf(alpha=pi.alpha, nn=nn, N=N, len=trunc.G, lseq=Ts, discount=discount)
      pi.prop        <- weights$pi.prop
      Vs             <- weights$Vs

    # Re-ordering & Slice Sampler
      index          <- order(pi.prop, decreasing=TRUE)
      pi.prop        <- pi.prop[index]
      Vs             <- Vs[index]
      mu             <- mu[,index, drop=FALSE]
      phi            <- phi[index]
      delta          <- delta[index]
      tau            <- tau[index]
      lmat           <- lmat[index]
      Qs             <- Qs[index]
      psi.inv        <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.numeric(levels(z))[z]
      if(!ind.slice) {
        ksi          <- pi.prop
        log.ksi      <- log(ksi)
      }
      u.slice        <- runif(N, 0, ksi[z])
      G              <- max(vapply(Ns, function(i) sum(u.slice[i] < ksi), integer(1L)))
      Gs             <- seq_len(G)
      log.slice.ind  <- vapply(Gs, function(g) slice.logs[1 + (u.slice < ksi[g])] - log.ksi[g], numeric(N))

    # Cluster Labels
      psi            <- 1/psi.inv
      Qgs            <- Qs[Gs]
      Q0             <- Qgs  > 0
      Q1             <- Qgs == 1
      if(G > 1)   {
        sigma        <- lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
        log.check    <- capture.output(log.probs <- try(vapply(Gs, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N)), silent=TRUE))
        if(inherits(log.probs, "try-error")) {
          log.probs  <- vapply(Gs, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) make.positive.definite(sigma[[g]]) else make.positive.definite(sqrt(sigma[[g]])), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
        }
        z.res        <- gumbel_max(probs=log.probs, log.like=TRUE, slice=TRUE)
        z            <- z.res$z
      } else      {
        z            <- rep(1, N)
      }
      nn             <- tabulate(z, nbins=trunc.G)
      nn0            <- nn   > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)
      dat.g          <- lapply(Gs, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data         <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta          <- base::matrix(0, nrow=N, ncol=0)
        eta.tmp      <- lapply(Gs, function(g) eta[z == g,, drop=FALSE])
        lmat[Gs]     <- lapply(Gs, base::matrix, 0, nrow=P, ncol=0)
      } else {
        eta.tmp      <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0, nrow=ifelse(Q0[g], 0, nn[g]), ncol=Qs[g]))
        EtE          <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[Gs]     <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Ps, function(j) .sim_load_s(Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]][,j],
                               EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nrow=P, byrow=TRUE) else
                               base::matrix(unlist(lapply(Ps, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nrow=P, byrow=FALSE))
      }

    # Uniquenesses
      psi.inv[,Gs]   <- vapply(Gs, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                               P=P, eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))

    # Means
      sum.data       <- vapply(dat.g, colSums, numeric(P))
      sum.eta        <- lapply(eta.tmp, colSums)
      mu[,Gs]        <- vapply(Gs, function(g) if(nn0[g]) .sim_mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                               sum.data=sum.data[,g], lmat=lmat[[g]], mu.zero=mu.zero) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Local Shrinkage
      load.2         <- lapply(lmat[Gs], .power2)
      phi[Gs]        <- lapply(Gs, function(g) if(nn0[g]) .sim_phi(Q=Qs[g], P=P, nu=nu, plus1=nuplus1,
                        tau=tau[[g]], load.2=load.2[[g]]) else .sim_phi_p(Q=Qs[g], P=P, nu=nu, plus1=nuplus1))

    # Global Shrinkage
      sum.terms      <- lapply(Gs, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gs) {
        Qg           <- Qs[g]
        Q1g          <- Q1[g]
        if(nn0[g]) {
          for(k in seq_len(Qg)) {
            delta[[g]][k]     <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P=P,
                                 Q=Qg, k=k, sum.term.kq=sum.terms[[g]][k:Qg]) else .sim_delta1(Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]],
                                 alpha.d1=ifelse(Q1g, alpha.d2, alpha.d1), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1])
            tau[[g]]          <- cumprod(delta[[g]])
          }
        } else {
          for(k in seq_len(Qg)) {
            delta[[g]][k]     <- if(k > 1) .sim_delta_p(alpha=alpha.d2, beta=beta.d2) else .sim_delta_p(alpha=ifelse(Q1g, alpha.d2, alpha.d1), beta=ifelse(Q1g, beta.d2, beta.d1))
            tau[[g]]          <- cumprod(delta[[g]])
          }
        }
      }

    # Adaptation
      if(all(adapt, iter > adapt.at)) {
        if(runif(1)   < ifelse(iter < burnin, 0.5, 1/exp(b0 + b1 * (iter - adapt.at)))) {
          colvec     <- lapply(nn.ind, function(g) (if(Q0[g]) colSums(abs(lmat[[g]]) < epsilon)/P else 0) >= prop)
          nonred     <- lapply(colvec, .which0)
          numred     <- lengths(colvec) - lengths(nonred)
          notred     <- numred == 0
          ng.ind     <- seq_along(nn.ind)
          Qs.old     <- Qs[nn0]
          Qs[nn0]    <- vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1 else Qs.old[h] - numred[h], numeric(1L))
          Q.big      <- Qs[nn0] > Q.star
          Q.bigs     <- any(Q.big)
          if(Q.bigs) {
            notred   <- notred & !Q.big
            Qs[nn0][Q.big]    <- Q.star
          }
          phi[nn0]   <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])],  rgamma(n=P, shape=nu + nuplus1, rate=nu)) else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0] <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])],     rgamma(n=1, shape=alpha.d2, rate=beta.d2)) else delta[[g]][nonred[[h]]])
          tau[nn0]   <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]  <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp       <- Qs[!nn0]
          Qpop       <- Qs[nn0]
          Fmax       <- max(Qpop)
          Qmax       <- ifelse(all(Q.big), Fmax, max(Qpop[!Q.big]))
          Qmaxold    <- max(Qs.old)
          if(Qmax     < max(Qemp, 0))  {
            Qs[Qmax   < Qs & !nn0] <- Qmax
            Qmaxseq  <- seq_len(Qmax)
            for(t in Ts[!nn0][Qemp  > Qmax]) {
              phi[[t]]        <- phi[[t]][,Qmaxseq,  drop=FALSE]
              delta[[t]]      <- delta[[t]][Qmaxseq]
              tau[[t]]        <- tau[[t]][Qmaxseq]
              lmat[[t]]       <- lmat[[t]][,Qmaxseq, drop=FALSE]
            }
          }
          if(all(sw["s.sw"], storage)) {
            eta.tmp  <- lapply(Gs,     function(g) if(nn0[g] && Qs[g] > Qs.old[which(nn.ind == g)]) cbind(eta.tmp[[g]], rnorm(nn[g])) else eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE])
          }
        }
      }

    # Alpha
      if(learn.alpha)      {
        if(discount  != 0) {
          MH.alpha   <- .sim_alpha_m(alpha=pi.alpha, discount=discount, alpha.shape=alpha.shape, alpha.rate=alpha.rate, N=N, G=G.non, zeta=zeta)
          pi.alpha   <- MH.alpha$alpha
          a.rate     <- MH.alpha$rate
        } else {
          pi.alpha   <- .sim_alpha_g(alpha=pi.alpha, shape=alpha.shape, rate=alpha.rate, G=G, N=N)
          a.rate     <- 1
        }
      }

    # Label Switching
      if(DP.lab.sw)   {
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
          } else        acc1  <- FALSE
        }
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
          } else        acc2  <- FALSE
        }
      }

      if(Q.bigs && !Q.large   && iter > burnin) {         warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large      <- TRUE
      }
      if(z.err  && !err.z) {                              warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z        <- TRUE
      }
      if(storage)    {
        if(verbose)      setTxtProgressBar(pb, iter)
        new.it       <-  which(iters == iter)
        if(sw["mu.sw"])  mu.store[,,new.it]            <- mu
        if(all(sw["s.sw"],
           any(Q0)))    {
          max.Q      <-  max(Qs)
          eta.tmp    <-  if(length(unique(Qs)) != 1)      lapply(Gs,       function(g) cbind(eta.tmp[[g]], base::matrix(0, nrow=nn[g], ncol=max.Q - Qs[g]))) else eta.tmp
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
        if(sw["pi.sw"])  pi.store[,new.it]             <- pi.prop
        if(learn.alpha)  alpha.store[new.it]           <- pi.alpha
        if(learn.d)      d.store[new.it]               <- discount
        if(MH.step)      a.rates[new.it]               <- a.rate
        if(learn.d)      d.rates[new.it]               <- d.rate
        if(DP.lab.sw)    lab.rate[,new.it]             <- c(acc1, acc2)
                         z.store[,new.it]              <- z
                         ll.store[new.it]              <- sum(z.res$log.like)
                         Q.store[,new.it]              <- Qs
                         G.store[new.it]               <- G.non
                         act.store[new.it]             <- G
      }
    }
    if(verbose)         close(pb)
    Gmax             <- seq_len(max(as.numeric(z.store)))
    Qmax             <- seq_len(max(Q.store))
    mu.store         <- if(sw["mu.sw"])  tryCatch(mu.store[,Gmax,, drop=FALSE],        error=function(e) mu.store)
    eta.store        <- if(sw["s.sw"])   tryCatch(eta.store[,Qmax,, drop=FALSE],       error=function(e) eta.store)
    load.store       <- if(sw["l.sw"])   tryCatch(load.store[,Qmax,Gmax,, drop=FALSE], error=function(e) load.store)
    psi.store        <- if(sw["psi.sw"]) tryCatch(psi.store[,Gmax,, drop=FALSE],       error=function(e) psi.store)
    pi.store         <- if(sw["pi.sw"])  tryCatch(pi.store[Gmax,, drop=FALSE],         error=function(e) pi.store)
    returns          <- list(mu        = if(sw["mu.sw"])  tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                             eta       = if(sw["s.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                             load      = if(sw["l.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                             psi       = if(sw["psi.sw"]) tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                             pi.prop   = if(sw["pi.sw"])  pi.store,
                             alpha     = if(learn.alpha)  alpha.store,
                             discount  = if(learn.d)      discount,
                             a.rate    = if(MH.step)      mean(a.rates) else a.rates,
                             d.rate    = if(learn.d)      mean(d.rates) else d.rates,
                             lab.rate  = if(DP.lab.sw)    setNames(rowmeans(lab.rate), c("Move1", "Move2")),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             Q.store   = tryCatch(Q.store[Gmax,, drop=FALSE],          error=function(e) Q.store),
                             G.store   = G.store,
                             act.store = act.store,
                             time      = init.time)
    attr(returns, "Q.big")    <- Q.large
    return(returns)
  }
