#########################################################################
### Gibbs Sampler for (Finite) Mixtures of Infinite Factor Analysers ####
#########################################################################

# Gibbs Sampler Function
  .gibbs_MIFA      <- function(Q, data, iters, N, P, G, sw, mu, mu.zero, uni.type, uni.prior,
                               col.mean, sigma.mu, burnin, thinning, verbose, nu1, nu2, cluster,
                               psi.alpha, psi.beta, adapt, truncated, start.AGS, stop.AGS, prop,
                               b0, b1, cluster.shrink, epsilon, equal.pro, forceQg, ...) {

  # Define & initialise variables
    start.time     <- proc.time()
    sq_mat         <- if(P   > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix         <- base::matrix
    total          <- max(iters)
    if(verbose)       pb    <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    AGS.burn       <- total/5L
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    Nseq           <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    uni            <- P == 1
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
    ll.store       <- integer(n.store)
    Q.store        <- matrix(0L, nrow=G, ncol=n.store)
    Q.large        <- Q.big <- Q.bigs <- FALSE
    err.z          <- z.err <- FALSE
    nu1.5          <- nu1 + 0.5
    P.5            <- P/2

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    if(all(mu.zero == 0)) {
      mu.zero      <- matrix(0L, nrow=1L, ncol=G)
      cluster$l.switch[1L]  <- FALSE
    }
    if(length(mu.zero)  == 1) {
      mu.zero      <- matrix(mu.zero, nrow=1L, ncol=G)
    }
    mu.prior       <- mu.sigma * mu.zero
    z              <- cluster$z
    nn             <- tabulate(z, nbins=G)
    nn0            <- nn > 0
    nn.ind         <- which(nn0)
    z.temp         <- factor(z, levels=Gseq)
    Q.star         <- Q
    Qs             <- rep(Q, G)
    Qs             <- if(forceQg) pmin(Qs, replace(nn, !nn0, Inf) - 1L) else Qs
    Q0             <- Qs  > 0
    Q1             <- Qs == 1
    pi.prop        <- cluster$pi.prop
    log.pis        <- log(pi.prop)
    pi.alpha       <- cluster$pi.alpha
    one.uni        <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv   <- switch(EXPR=uni.type,  unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                             constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip    <- switch(EXPR=uni.prior, unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni))       {
      uni.shape    <- switch(EXPR=uni.type,  constrained=N/2 + psi.alpha, single=(N * P)/2 + psi.alpha)
      V            <- switch(EXPR=uni.type,  constrained=P, single=1L)
    }
    if(uni.prior   == "isotropic")   {
      psi.beta     <- matrix(vapply(Gseq, function(g) psi.beta[which.max(.ndeci(psi.beta[,g])),g], numeric(1L)), nrow=1, ncol=G)
    } else if(length(psi.beta) == 1) {
      psi.beta     <- matrix(psi.beta, nrow=1L, ncol=G)
    }
    alpha.d1       <- cluster$alpha.d1
    alpha.d2       <- cluster$alpha.d2
    ad1.x          <- length(unique(alpha.d1)) == 1L
    ad2.x          <- length(unique(alpha.d2)) == 1L
    beta.d1        <- cluster$beta.d1
    beta.d2        <- cluster$beta.d2
    bd1.x          <- length(unique(beta.d1))  == 1L
    bd2.x          <- length(unique(beta.d2))  == 1L
    rho1           <- cluster$rho1
    rho2           <- cluster$rho2
    r1.x           <- length(unique(rho1))     == 1L
    r2.x           <- length(unique(rho2))     == 1L
    mu0g           <- cluster$l.switch[1L]
    psi0g          <- cluster$l.switch[2L]
    delta0g        <- cluster$l.switch[3L]
    label.switch   <- any(cluster$l.switch)
    Qmax           <- ifelse(forceQg, max(Qs), Q)
    Qmaxseq        <- seq_len(Qmax)
    eta            <- .sim_eta_p(N=N, Q=Qmax)
    phi            <- if(forceQg) lapply(Gseq, function(g) .sim_phi_p(Q=Qs[g], P=P, nu1=nu1, nu2=nu2)) else replicate(G, .sim_phi_p(Q=Q, P=P, nu1=nu1, nu2=nu2), simplify=FALSE)
    if(isTRUE(truncated))       {
      .sim_deltak  <- .sim_deltaKT
      .rdelta      <- rltrgamma
      delta        <- lapply(Gseq, function(g) c(if(!forceQg | Qs[g] > 0) .sim_delta_p(alpha=alpha.d1[g], beta=beta.d1[g]), .sim_deltaPT(Q=Qs[g], alpha=alpha.d2[g], beta=beta.d2[g])))
      .sim_delta_p <- .sim_deltaPT
    } else          {
      .rdelta      <- stats::rgamma
      delta        <- lapply(Gseq, function(g) c(if(!forceQg | Qs[g] > 0) .sim_delta_p(alpha=alpha.d1[g], beta=beta.d1[g]), .sim_delta_p(Q=Qs[g], alpha=alpha.d2[g], beta=beta.d2[g])))
    }
    tau            <- lapply(delta, cumprod)
    if(cluster.shrink)  {
      sig.store    <- matrix(0L, nrow=G, ncol=n.store)
      MGPsig       <- .sim_sigma_p(G=G, rho1=rho1, rho2=rho2)
    } else MGPsig  <- rep(1L, G)
    lmat           <- lapply(Gseq, function(g) matrix(vapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g])), nrow=P, byrow=TRUE))
    if(isTRUE(one.uni)) {
      psi.beta     <- psi.beta[,1L]
      psi.inv      <- matrix(.sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), nrow=P, ncol=G)
    } else psi.inv <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    psi.inv        <- if(uni)     t(psi.inv)    else psi.inv
    if(isTRUE(one.uni)) {
      psi.inv[]    <- 1/switch(EXPR=uni.type, constrained=colVars(data, center=col.mean, refine=FALSE, useNames=FALSE), max(colVars(data, center=col.mean, refine=FALSE, useNames=FALSE)))
    } else   {
      tmp.psi      <- (nn[nn0] - 1L)/pmax(rowsum(data^2, z) - rowsum(data, z)^2/nn[nn0], 0L)
      tmp.psi      <- switch(EXPR=uni.type, unconstrained=t(tmp.psi), matrix(Rfast::rowMaxs(tmp.psi, value=TRUE), nrow=P, ncol=G, byrow=TRUE))
      psi.inv[,nn   > 1]    <- tmp.psi[!is.nan(tmp.psi)]
      rm(tmp.psi)
    }
    max.p          <- (psi.alpha  - 1)/psi.beta
    inf.ind        <- psi.inv > max(max.p)
    psi.inv[inf.ind]        <- matrix(max.p, nrow=P, ncol=G)[inf.ind]
    rm(max.p, inf.ind)
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total))   {
      if(verbose   && iter   < burnin) utils::setTxtProgressBar(pb, iter)
      storage      <- is.element(iter, iters)

      # Adaptation
      if(adapt     && all(iter >= start.AGS, iter < stop.AGS))    {
        if(stats::runif(1) < ifelse(iter < AGS.burn, 0.5, exp(-b0 - b1 * (iter - start.AGS))))  {
          colvec   <- lapply(nn.ind, function(g) if(Q0[g]) (colSums2(abs(lmat[[g]])      < epsilon,
                                                                     useNames=FALSE)/P) >= prop else stats::runif(1) <= prop)
          nonred   <- lapply(colvec, .which0)
          numred   <- lengths(colvec) - lengths(nonred)
          notred   <- numred == 0
          ng.ind   <- seq_along(nn.ind)
          Qs.old   <- Qs[nn0]
          Qs[nn0]  <- pmax.int(0L, vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1L      else Qs.old[h] - numred[h], numeric(1L)))
          star.Q   <- if(forceQg) pmin(Q.star, nn[nn0] - 1L) else Q.star
          Q.big    <- Qs[nn0] > star.Q
          if((Q.bigs        <-  any(Q.big))) {
            notred <- notred & !Q.big
            Qs[nn0][Q.big]  <- if(forceQg) star.Q[Q.big]     else Q.star
            if(forceQg)      {
              for(qb in which(Q.big))   {
               nonred[[qb]] <- seq_len(star.Q[qb])
              }
            }
          }
          phi[nn0]          <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])],  .rgamma0(n=P, shape=nu1,         rate=nu2))        else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0]        <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])],     .rdelta(n=1,  shape=alpha.d2[g], rate=beta.d2[g])) else delta[[g]][nonred[[h]]])
          tau[nn0]          <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]         <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]] * MGPsig[g]))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp     <- Qs[!nn0]
          Qmax     <- max(Qs[nn0])
          Qmaxseq  <- seq_len(Qmax)
          if(any(!nn0)    && max(Qemp) != Qmax)  {
            for(g  in Gseq[!nn0][Qemp  != Qmax]) {
              Qg   <- Qs[g]
              if(Qg > Qmax)  {
                phi[[g]]    <- phi[[g]][,Qmaxseq,  drop=FALSE]
                delta[[g]]  <- delta[[g]][Qmaxseq]
                tau[[g]]    <- tau[[g]][Qmaxseq]
                lmat[[g]]   <- lmat[[g]][,Qmaxseq, drop=FALSE]
              } else {
                while(Qg    != Qmax)   {
                 phi[[g]]   <- cbind(phi[[g]],  .rgamma0(n=P, shape=nu1,         rate=nu2))
                 delta[[g]] <- c(delta[[g]],    .rdelta(n=1,  shape=alpha.d2[g], rate=beta.d2[g]))
                 tau[[g]]   <- cumprod(delta[[g]])
                 Qg         <- Qg + 1L
                 lmat[[g]]  <- cbind(lmat[[g]], stats::rnorm(n=P, mean=0, sd=1/sqrt(phi[[g]][,Qg] * tau[[g]][Qg] * MGPsig[g])))
                }
              }
            }
            Qs[Qmax != Qs    & !nn0] <- Qmax
          }
          Q0        <- Qs  > 0
          Q1        <- Qs == 1
        }
      }

    # Mixing Proportions
      pi.prop      <- if(equal.pro) pi.prop else rDirichlet(G=G, alpha=pi.alpha, nn=nn)

    # Scores & Loadings
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2L, mu[,g], FUN="-", check.margin=FALSE))
      n.eta        <- nn
      n0q0         <- nn0  & Q0
      q0ng         <- (!Q0 | Q1) & n.eta > 0
      if(!any(Q0))    {
        eta        <- .empty_mat(nr=N)
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,,  drop=FALSE])
        lmat       <- replicate(G, .empty_mat(nr=P))
      } else {
        eta.tmp    <- lapply(Gseq, function(g) if(n0q0[g]) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0L, nrow=ifelse(Q0[g], 0L, nn[g]), ncol=Qs[g]))
        EtE        <- lapply(Gseq, function(g) if(n0q0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- lapply(Gseq, function(g) matrix(if(n0q0[g]) vapply(Pseq, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j], Q1=Q1[g],
                             EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g]))     else
                             vapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], sigma=MGPsig[g]), numeric(Qs[g])), nrow=P, byrow=TRUE))
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat      <- lapply(Gseq, function(g) { S <- c.data[[g]] - if(Q0[g]) tcrossprod(eta.tmp[[g]], lmat[[g]]) else 0L; S^2 } )
        psi.inv[,] <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,] <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta[,g], lmat=lmat[[g]],
                      P=P, Q0=Q0[g], eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
      }

    # Means
      sum.data     <- vapply(dat.g, colSums2, useNames=FALSE, numeric(P))
      sum.data     <- if(uni) t(sum.data) else sum.data
      sum.eta      <- lapply(eta.tmp, colSums2, useNames=FALSE)
      mu[,]        <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.prior=mu.prior[,g], sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                             sum.data=sum.data[,g], lmat=lmat[[g]], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero[,g]), numeric(P))

    # Shrinkage
      if(any(Q0))     {
        load.2     <- lapply(lmat, "^", 2)
        phi        <- lapply(Gseq, function(g) if(n0q0[g]) .sim_phi(Q=Qs[g], P=P, nu1.5=nu1.5, nu2=nu2, tau=tau[[g]],
                      load.2=load.2[[g]], sigma=MGPsig[g]) else .sim_phi_p(Q=Qs[g], P=P, nu1=nu1, nu2=nu2))
        sum.terms  <- lapply(Gseq, function(g) if(n0q0[g]) colSums2(phi[[g]] * load.2[[g]], useNames=FALSE))
        for(g in Gseq)  {
          Qg       <- Qs[g]
          Q1g      <- Q1[g]
          if(n0q0[g])   {
            for(k in seq_len(Qg)) {
              delta[[g]][k] <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2[g], beta.d2=beta.d2[g], delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P.5=P.5, Q=Qg,
                               k=k, sum.term.kq=sum.terms[[g]][k:Qg], sigma=MGPsig[g]) else .sim_delta1(Q=Qg, P.5=P.5, tau=tau[[g]], sum.term=sum.terms[[g]],
                               alpha.d1=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), beta.d1=ifelse(Q1g, beta.d2[g], beta.d1[g]), delta.1=delta[[g]][1L], sigma=MGPsig[g])
              tau[[g]]      <- cumprod(delta[[g]])
            }
          } else {
            if(Q0[g])   {
              delta[[g]]    <- c(stats::rgamma(n=1, shape=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), rate=ifelse(Q1g, beta.d2[g], beta.d1[g])),
                                 .sim_delta_p(Q=Qg, alpha=alpha.d2[g], beta=beta.d2[g]))
              tau[[g]]      <- cumprod(delta[[g]])
            }
          }
        }
        if(cluster.shrink)   {
          nGq0              <- sum(n0q0)
          MGPsig[n0q0]      <- .sim_sigma(G=nGq0, P.5=P.5, Qs=Qs[n0q0], rho1=rho1[n0q0], rho2=rho2[n0q0], sum.terms=sum.terms[n0q0], tau=tau[n0q0])
          MGPsig[!n0q0]     <- .sim_sigma_p(G=G - nGq0, rho1=rho1[n0q0], rho2=rho2[n0q0])
        }
      }

    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- if(uni) lapply(Gseq, function(g) as.matrix(psi[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0L)) else lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      log.pis      <- if(equal.pro) log.pis else log(pi.prop)
      if(uni) {
        log.probs  <- vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log.pis[g], numeric(N))
      } else  {
        log.probs  <- try(vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sq_mat(sigma[[g]]), log=TRUE, isChol=!Q) + log.pis[g], numeric(N)), silent=TRUE)
        if(z.err   <- inherits(log.probs, "try-error")) {
         log.probs <- vapply(Gseq, function(g, Q=Q0[g]) { sigma <- if(Q) is.posi_def(sigma[[g]], make=TRUE)$X.new else sq_mat(sigma[[g]]); dmvn(data, mu[,g], sigma, log=TRUE, isChol=!Q) + log.pis[g] }, numeric(N))
        }
      }
      z            <- gumbel_max(probs=log.probs)

    # Label Switching
      if(label.switch)   {
        sw.lab     <- .lab_switch(z.new=z, z.old=z.temp)
        z.perm     <- sw.lab$z.perm
        left       <- as.integer(unname(z.perm))
        right      <- as.integer(names(z.perm))
        if(!identical(left, right))   {
          z        <- sw.lab$z
          if(length(unique(Qs)) != 1) {
            Qs[left]        <- Qs[right]
          }
          mu[,left]         <- mu[,right,       drop=FALSE]
          lmat[left]        <- lmat[right]
          delta[left]       <- delta[right]
          phi[left]         <- phi[right]
          tau[left]         <- tau[right]
          psi.inv[,left]    <- psi.inv[,right,  drop=FALSE]
          pi.prop[left]     <- pi.prop[right]
          nn[left]          <- nn[right]
          Q0[left]          <- Q0[right]
          Q1[left]          <- Q1[right]
          if(mu0g)           {
            mu.zero[,left]  <- mu.zero[,right,  drop=FALSE]
            mu.prior[,left] <- mu.prior[,right, drop=FALSE]
          }
          if(psi0g)          {
            psi.beta[,left] <- psi.beta[,right, drop=FALSE]
          }
          if(all(delta0g,
                 !ad1.x))    {
            alpha.d1[left]  <- alpha.d1[right]
          }
          if(all(delta0g,
                 !ad2.x))    {
            alpha.d2[left]  <- alpha.d2[right]
          }
          if(all(delta0g,
                 !bd1.x))    {
            beta.d1[left]   <- beta.d1[right]
          }
          if(all(delta0g,
                 !bd2.x))    {
            beta.d2[left]   <- beta.d2[right]
          }
          if(cluster.shrink) {
             MGPsig[left]   <- MGPsig[right]
             if(delta0g)     {
               if(!r1.x)     {
                 rho1[left] <- rho1[right]
               }
               if(!r2.x)     {
                 rho2[left] <- rho2[right]
               }
             }
          }
        }
      }
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn  > 0
      nn.ind       <- which(nn0)

      if(Q.bigs && !Q.large && iter > burnin) {        cat("\n"); warning(paste0("\nQ has exceeded initial number of loadings columns", ifelse(forceQg, " (or exceeded the number of observations in one or more clusters)", ""), " since burnin:\nconsider increasing 'range.Q' from ", Q.star, ifelse(forceQg, " or setting 'forceQg' to FALSE\n", "\n")), call.=FALSE, immediate.=TRUE)
        Q.large    <- TRUE
      }
      if(z.err  && !err.z) {                           cat("\n"); warning("\nAlgorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices\n", call.=FALSE, immediate.=TRUE)
        err.z      <- TRUE
      }
      if(storage)   {
        if(verbose)   utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])        mu.store[,,new.it]   <-   mu
        if(all(sw["s.sw"],
           any(Q0)))     {
          eta.tmp  <- if(length(unique(Qs)) != 1)        lapply(Gseq,       function(g) cbind(eta.tmp[[g]], base::matrix(0L, nrow=n.eta[g], ncol=Qmax - ncol(eta.tmp[[g]])))) else eta.tmp
          if(any(q0ng))  {
            eta.tmp[q0ng]   <-                           lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                       eta.store[,Qmaxseq,new.it]   <-   do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])   {
          for(g in Gseq) {
            Qseqg  <- seq_len(Qs[g])
            if(Q0[g]) load.store[,Qseqg,g,new.it]   <-   lmat[[g]]
          }
        }
        if(sw["psi.sw"])      psi.store[,,new.it]   <-   1/psi.inv
        if(sw["pi.sw"])         pi.store[,new.it]   <-   pi.prop
        if(cluster.shrink)     sig.store[,new.it]   <-   MGPsig
                                 z.store[new.it,]   <-   as.integer(z)
                                 ll.store[new.it]   <-   sum(rowLogSumExps(log.probs, useNames=FALSE))
                                 Q.store[,new.it]   <-   as.integer(Qs)
      }
    }
    if(verbose)       close(pb)
    Qmax           <- seq_len(max(Q.store))
    eta.store      <- if(sw["s.sw"])  tryCatch(eta.store[,Qmax,, drop=FALSE],   error=function(e) eta.store)
    load.store     <- if(sw["l.sw"])  tryCatch(load.store[,Qmax,,, drop=FALSE], error=function(e) load.store)
    returns        <- list(mu       = if(sw["mu.sw"])    tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                           eta      = if(sw["s.sw"])     tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                           load     = if(sw["l.sw"])     tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                           psi      = if(sw["psi.sw"])   tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                           pi.prop  = if(sw["pi.sw"])    pi.store,
                           sigma    = if(cluster.shrink) sig.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           Q.store  = Q.store,
                           time     = init.time)
    attr(returns, "Q.big")  <- Q.large
      return(returns)
  }
