###########################################################################
### Gibbs Sampler for Overfitted Mixtures of Infinite Factor Analysers ####
###########################################################################

# Gibbs Sampler Function
  .gibbs_OMIFA       <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu, uni.type, uni.prior, truncated, burnin,
                                 thinning, adapt, psi.alpha, psi.beta, verbose, alpha.d1, alpha.d2, sw, cluster,
                                 nu1, nu2, rho1, rho2, b0, b1, mu, prop, beta.d1, beta.d2, start.AGS, stop.AGS,
                                 epsilon, learn.alpha, a.hyper, zeta, tune.zeta, cluster.shrink, forceQg, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    sq_mat           <- if(P   > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix           <- base::matrix
    total            <- max(iters)
    if(verbose)         pb    <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    AGS.burn         <- total/5L
    Gseq             <- seq_len(G)
    Pseq             <- seq_len(P)
    Nseq             <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
    uni              <- P == 1
    if(sw["mu.sw"])  {
      mu.store       <- array(0L,  dim=c(P, G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store      <- array(0L,  dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store     <- array(0L,  dim=c(P, Q, G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store      <- array(0L,  dim=c(P, G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store       <- matrix(0L, nrow=G, ncol=n.store)
    }
    z.store          <- matrix(0L, nrow=n.store, ncol=N)
    ll.store         <- vector("integer", n.store)
    Q.store          <- matrix(0L, nrow=G, ncol=n.store)
    Q.large          <- Q.big <- Q.bigs <- FALSE
    err.z            <- z.err <- FALSE
    G.store          <- vector("integer", n.store)
    nu1.5            <- nu1 + 0.5
    P.5              <- P/2

    mu.sigma         <- 1/sigma.mu
    mu.prior         <- mu.sigma * mu.zero
    sig.mu.sqrt      <- sqrt(sigma.mu)
    z                <- cluster$z
    nn               <- tabulate(z, nbins=G)
    nn0              <- nn  > 0
    nn.ind           <- which(nn0)
    G.non            <- length(nn.ind)
    Q.star           <- Q
    Qs               <- rep(Q, G)
    Qs               <- if(forceQg) pmin(Qs, replace(nn, !nn0, Inf) - 1L) else Qs
    Q0               <- Qs  > 0
    Q1               <- Qs == 1
    pi.alpha         <- cluster$pi.alpha
    if(learn.alpha)   {
      alpha.store    <- ll.store
      alpha.shape    <- a.hyper[1L]
      alpha.rate     <- a.hyper[2L]
      a.rates        <- vector("integer", total)
    } else a.rates   <- 1L
    avgzeta          <- zeta
    heat             <- tune.zeta$heat
    lambda           <- tune.zeta$lambda
    target           <- tune.zeta$target
    zeta.tune        <- tune.zeta$do
    startz           <- tune.zeta$start.zeta
    stopz            <- tune.zeta$stop.zeta
    one.uni          <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv     <- switch(EXPR=uni.type,  unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                               constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip      <- switch(EXPR=uni.prior, unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni))   {
      uni.shape      <- switch(EXPR=uni.type,  constrained=N/2 + psi.alpha, single=(N * P)/2 + psi.alpha)
      V              <- switch(EXPR=uni.type,  constrained=P, single=1L)
    }
    psi.beta         <- switch(EXPR=uni.prior, isotropic=psi.beta[which.max(.ndeci(psi.beta))], psi.beta)
    pi.prop          <- c(cluster$pi.prop, vector("numeric", G - length(cluster$pi.prop)))
    mu               <- cbind(mu, vapply(seq_len(G - length(cluster$pi.prop)), function(g) .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P)))
    Qmax             <- ifelse(forceQg, max(Qs), Q)
    Qmaxseq          <- seq_len(Qmax)
    eta              <- .sim_eta_p(N=N, Q=Qmax)
    phi              <- if(forceQg) lapply(Gseq, function(g) .sim_phi_p(Q=Qs[g], P=P, nu1=nu1, nu2=nu2)) else replicate(G, .sim_phi_p(Q=Q, P=P, nu1=nu1, nu2=nu2), simplify=FALSE)
    if(isTRUE(truncated))   {
      .sim_deltak    <- .sim_deltaKT
      .rdelta        <- rltrgamma
      delta          <- if(forceQg) lapply(Gseq, function(g) c(if(Qs[g] > 0) .sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_deltaPT(Q=Qs[g], alpha=alpha.d2, beta=beta.d2))) else replicate(G, list(c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_deltaPT(Q=Q, alpha=alpha.d2, beta=beta.d2))))
      .sim_delta_p   <- .sim_deltaPT
    } else            {
      .rdelta        <- stats::rgamma
      delta          <- if(forceQg) lapply(Gseq, function(g) c(if(Qs[g] > 0) .sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Qs[g], alpha=alpha.d2, beta=beta.d2))) else replicate(G, list(c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2))))
    }
    if(cluster.shrink)    {
      sig.store      <- matrix(0L, nrow=G, ncol=n.store)
      MGPsig         <- .sim_sigma_p(G=G, rho1=rho1, rho2=rho2)
    } else MGPsig    <- rep(1L, G)
    tau              <- lapply(delta, cumprod)
    lmat             <- lapply(Gseq, function(g) matrix(vapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g])), nrow=P, byrow=TRUE))
    if(isTRUE(one.uni)) {
      psi.inv        <- matrix(.sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), nrow=P, ncol=G)
    } else psi.inv   <- replicate(G, .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), simplify="array")
    psi.inv          <- if(uni) t(psi.inv) else psi.inv
    if(isTRUE(one.uni)) {
      psi.inv[]      <- 1/switch(EXPR=uni.type, constrained=colVars(data), max(colVars(data)))
    } else {
      tmp.psi        <- (nn[nn0] - 1L)/pmax(rowsum(data^2, z) - rowsum(data, z)^2/nn[nn0], 0L)
      tmp.psi        <- switch(EXPR=uni.type, unconstrained=t(tmp.psi), matrix(Rfast::rowMaxs(tmp.psi, value=TRUE), nrow=P, ncol=G, byrow=TRUE))
      psi.inv[,nn     > 1]    <- tmp.psi[!is.nan(tmp.psi)]
      rm(tmp.psi)
    }
    max.p            <- (psi.alpha  - 1)/psi.beta
    inf.ind          <- psi.inv > max(max.p)
    psi.inv[inf.ind] <- matrix(max.p, nrow=P, ncol=G)[inf.ind]
    rm(max.p, inf.ind)
    init.time        <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total))     {
      if(verbose     && iter   < burnin) utils::setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

      # Adaptation
      if(adapt       && all(iter >= start.AGS, iter < stop.AGS))  {
        if(stats::runif(1) < ifelse(iter < AGS.burn, 0.5, exp(-b0 - b1 * (iter - start.AGS))))    {
          colvec     <- lapply(nn.ind, function(g) (if(Q0[g]) colSums(abs(lmat[[g]])    < epsilon)/P else stats::runif(1)) >= prop)
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
            for(g in Gseq[!nn0][Qemp    != Qmax]) {
              Qg     <- Qs[g]
              if(Qg   > Qmax)  {
                phi[[g]]      <- phi[[g]][,Qmaxseq,  drop=FALSE]
                delta[[g]]    <- delta[[g]][Qmaxseq]
                tau[[g]]      <- tau[[g]][Qmaxseq]
                lmat[[g]]     <- lmat[[g]][,Qmaxseq, drop=FALSE]
              } else  {
                while(Qg   != Qmax)   {
                 phi[[g]]     <- cbind(phi[[g]],     .rgamma0(n=P, shape=nu1,      rate=nu2))
                 delta[[g]]   <- c(delta[[g]],       .rdelta(n=1,  shape=alpha.d2, rate=beta.d2))
                 tau[[g]]     <- cumprod(delta[[g]])
                 Qg  <- Qg  + 1L
                 lmat[[g]]    <- cbind(lmat[[g]],    stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[[g]][,Qg] * tau[[g]][Qg] * MGPsig[g])))
                }
              }
            }
            Qs[Qmax  != Qs  & !nn0]  <- Qmax
          }
          Q0         <- Qs  > 0
          Q1         <- Qs == 1
        }
      }

    # Mixing Proportions & Re-ordering
      pi.prop        <- rDirichlet(G=G, alpha=pi.alpha, nn=nn)
      index          <- order(nn, decreasing=TRUE)
      pi.prop        <- pi.prop[index]
      mu             <- mu[,index, drop=FALSE]
      phi            <- phi[index]
      delta          <- delta[index]
      tau            <- tau[index]
      MGPsig         <- MGPsig[index]
      lmat           <- lmat[index]
      Qs             <- Qs[index]
      Q0             <- Q0[index]
      Q1             <- Q1[index]
      psi.inv        <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.integer(levels(z))[z]
      nn             <- nn[index]
      nn0            <- nn0[index]

    # Scores & Loadings
      dat.g          <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])
      c.data         <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2L, mu[,g], FUN="-", check.margin=FALSE))
      n.eta          <- nn
      n0q0           <- nn0  & Q0
      q0ng           <- (!Q0 | Q1) & n.eta > 0
      if(!any(Q0))    {
        eta          <- .empty_mat(nr=N)
        eta.tmp      <- lapply(Gseq, function(g) eta[z == g,,  drop=FALSE])
        lmat         <- replicate(G, .empty_mat(nr=P))
      } else {
        eta.tmp      <- lapply(Gseq, function(g) if(n0q0[g]) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0L, nrow=ifelse(Q0[g], 0L, nn[g]), ncol=Qs[g]))
        EtE          <- lapply(Gseq, function(g) if(n0q0[g]) crossprod(eta.tmp[[g]]))
        lmat         <- lapply(Gseq, function(g) matrix(if(n0q0[g]) vapply(Pseq, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j], Q1=Q1[g],
                               EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g]))     else
                               vapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g])), nrow=P, byrow=TRUE))
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat        <- lapply(Gseq, function(g) { S <- c.data[[g]] - if(Q0[g]) tcrossprod(eta.tmp[[g]], lmat[[g]]) else 0L; S^2 } )
        psi.inv[,]   <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,]   <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                        P=P, Q0=Q0[g], eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      }

    # Means
      sum.data       <- vapply(dat.g, colSums2, numeric(P))
      sum.data       <- if(uni) t(sum.data) else sum.data
      sum.eta        <- lapply(eta.tmp, colSums2)
      mu[,]          <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                               sum.data=sum.data[,g], lmat=lmat[[g]], mu.prior=mu.prior) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Shrinkage
      if(any(Q0))    {
        load.2       <- lapply(lmat, "^", 2)
        phi          <- lapply(Gseq, function(g) if(n0q0[g]) .sim_phi(Q=Qs[g], P=P, nu1.5=nu1.5, nu2=nu2, tau=tau[[g]],
                        load.2=load.2[[g]], sigma=MGPsig[g]) else .sim_phi_p(Q=Qs[g], P=P, nu1=nu1, nu2=nu2))
        sum.terms    <- lapply(Gseq, function(g) if(n0q0[g]) colSums2(phi[[g]] * load.2[[g]]))
        for(g in Gseq) {
          Qg         <- Qs[g]
          Q1g        <- Q1[g]
          if(n0q0[g])  {
            for(k in seq_len(Qg)) {
              delta[[g]][k]   <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P.5=P.5, Q=Qg,
                                 k=k, sum.term.kq=sum.terms[[g]][k:Qg], sigma=MGPsig[g]) else .sim_delta1(Q=Qg, P.5=P.5, tau=tau[[g]], sum.term=sum.terms[[g]],
                                 alpha.d1=ifelse(Q1g, alpha.d2, alpha.d1), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1L], sigma=MGPsig[g])
              tau[[g]]        <- cumprod(delta[[g]])
            }
          } else {
            if(Q0[g])  {
              delta[[g]]      <- c(stats::rgamma(n=1, shape=ifelse(Q1g, alpha.d2, alpha.d1), rate=ifelse(Q1g, beta.d2, beta.d1)),
                                   .sim_delta_p(Q=Qg, alpha=alpha.d2, beta=beta.d2))
              tau[[g]]        <- cumprod(delta[[g]])
            }
          }
        }
        if(cluster.shrink)     {
          nGq0                <- sum(n0q0)
          MGPsig[n0q0]        <- .sim_sigma(G=nGq0, P.5=P.5, Qs=Qs[n0q0], rho1=rho1, rho2=rho2, sum.terms=sum.terms[n0q0], tau=tau[n0q0])
          MGPsig[!n0q0]       <- .sim_sigma_p(G=G - nGq0, rho1=rho1, rho2=rho2)
        }
      }

    # Cluster Labels
      psi            <- 1/psi.inv
      sigma          <- if(uni) lapply(Gseq, function(g) as.matrix(psi[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0L)) else lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      log.pis        <- log(pi.prop)
      if(uni) {
        log.probs    <- vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log.pis[g], numeric(N))
      } else  {
        log.probs    <- try(vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sq_mat(sigma[[g]]), log=TRUE, isChol=!Q) + log.pis[g], numeric(N)), silent=TRUE)
        if(z.err     <- inherits(log.probs, "try-error")) {
          log.probs  <- vapply(Gseq, function(g, Q=Q0[g]) { sigma <- if(Q) is.posi_def(sigma[[g]], make=TRUE)$X.new    else sq_mat(sigma[[g]]); dmvn(data, mu[,g], sigma, log=TRUE, isChol=!Q) + log.pis[g] }, numeric(N))
        }
      }
      z              <- gumbel_max(probs=log.probs)
      nn             <- tabulate(z, nbins=G)
      nn0            <- nn  > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)

    # Alpha
      if(learn.alpha)      {
        MH.alpha     <- .sim_alpha_o(alpha=pi.alpha, zeta=zeta, G=G, N=N, nn=nn[nn0], shape=alpha.shape, rate=alpha.rate)
        pi.alpha     <- MH.alpha$alpha
        a.rates[iter]         <- MH.alpha$rate
        if(isTRUE(zeta.tune))   {
          if(iter    >= startz &&
             iter     < stopz)  {
            zeta     <- .tune_zeta(zeta=zeta, time=iter, l.rate=MH.alpha$l.prob, heat=heat, target=target, lambda=lambda)
          }
          avgzeta    <- c(avgzeta, zeta)
        }
      }

      if(Q.bigs && !Q.large   && iter > burnin) {         cat("\n"); cat("\n"); warning(paste0("\nQ has exceeded initial number of loadings columns", ifelse(forceQg, " (or exceeded the number of observations in one or more clusters)", ""), " since burnin:\nconsider increasing 'range.Q' from ", Q.star, ifelse(forceQg, " or setting 'forceQg' to FALSE\n", "\n")), call.=FALSE)
        Q.large      <- TRUE
      }
      if(z.err  && !err.z) {                              cat("\n"); warning("\nAlgorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices\n", call.=FALSE)
        err.z        <- TRUE
      }
      if(storage)    {
        if(verbose)      utils::setTxtProgressBar(pb, iter)
        new.it       <-  which(iters == iter)
        if(sw["mu.sw"])           mu.store[,,new.it]   <-    mu
        if(all(sw["s.sw"],
           any(Q0)))     {
          eta.tmp    <-  if(length(unique(Qs)) != 1)         lapply(Gseq,       function(g) cbind(eta.tmp[[g]], base::matrix(0L, nrow=n.eta[g], ncol=Qmax - Qs[g]))) else eta.tmp
          if(any(q0ng))  {
            eta.tmp[q0ng]     <-                             lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                          eta.store[,Qmaxseq,new.it]   <-    do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])   {
          for(g in Gseq) {
            Qseqg    <-  seq_len(Qs[g])
            if(Q0[g])    load.store[,Qseqg,g,new.it]   <-    lmat[[g]]
          }
        }
        if(sw["psi.sw"])         psi.store[,,new.it]   <-    1/psi.inv
        if(sw["pi.sw"])            pi.store[,new.it]   <-    pi.prop
        if(learn.alpha)          alpha.store[new.it]   <-    pi.alpha
        if(cluster.shrink)        sig.store[,new.it]   <-    MGPsig
                                    z.store[new.it,]   <-    as.integer(z)
                                    ll.store[new.it]   <-    sum(rowLogSumExps(log.probs))
                                    Q.store[,new.it]   <-    as.integer(Qs)
                                     G.store[new.it]   <-    as.integer(G.non)
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
    returns          <- list(mu        = if(sw["mu.sw"])     tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                             eta       = if(sw["s.sw"])      tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                             load      = if(sw["l.sw"])      tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                             psi       = if(sw["psi.sw"])    tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                             pi.prop   = if(sw["pi.sw"])     pi.store,
                             sigma     = if(cluster.shrink)  sig.store,
                             alpha     = if(learn.alpha)     alpha.store,
                             a.rate    = ifelse(learn.alpha, mean(a.rates), a.rates),
                             z.store   = z.store,
                             ll.store  = ll.store,
                             G.store   = G.store,
                             Q.store   = tryCatch(Q.store[Gmax,, drop=FALSE],          error=function(e) Q.store),
                             avg.zeta  = if(learn.alpha)     ifelse(zeta.tune, mean(avgzeta), zeta),
                             time      = init.time)
    attr(returns, "Q.big")    <- Q.large
      return(returns)
  }
