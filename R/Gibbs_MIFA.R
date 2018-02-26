#########################################################################
### Gibbs Sampler for (Finite) Mixtures of Infinite Factor Analysers ####
#########################################################################

# Gibbs Sampler Function
  .gibbs_MIFA      <- function(Q, data, iters, N, P, G, sw, mu, mu.zero, uni.type,
                               uni.prior, sigma.mu, burnin, thinning, verbose, nu,
                               cluster, psi.alpha, psi.beta, adapt, adaptat, prop,
                               b0, b1, beta.d1, beta.d2, epsilon, nuplus1, equal.pro, ...) {

  # Define & initialise variables
    start.time     <- proc.time()
    sq_mat         <- if(P   > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix         <- base::matrix
    total          <- max(iters)
    if(verbose)       pb    <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
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
    ll.store       <- vector("integer", n.store)
    Q.star         <- Q
    Qs             <- rep(Q, G)
    Q.store        <- matrix(0L, nrow=G, ncol=n.store)
    Q.large        <- Q.big <- Q.bigs <- FALSE
    err.z          <- z.err <- FALSE

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    if(all(mu.zero == 0)) {
      mu.zero      <- matrix(0L, nrow=1, ncol=G)
      cluster$l.switch[1]   <- FALSE
    }
    if(length(mu.zero)  == 1) {
      mu.zero      <- matrix(mu.zero, nrow=1, ncol=G)
    }
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    nn0            <- nn > 0
    pi.prop        <- cluster$pi.prop
    log.pis        <- log(pi.prop)
    pi.alpha       <- cluster$pi.alpha
    one.uni        <- is.element(uni.type, c("constrained", "single"))
    .sim_psi_inv   <- switch(uni.type,   unconstrained=.sim_psi_uu,   isotropic=.sim_psi_uc,
                                         constrained=.sim_psi_cu,     single=.sim_psi_cc)
    .sim_psi_ip    <- switch(uni.prior,  unconstrained=.sim_psi_ipu,  isotropic=.sim_psi_ipc)
    if(isTRUE(one.uni))       {
      uni.shape    <- switch(uni.type,   constrained=N/2 + psi.alpha, single=(N * P)/2 + psi.alpha)
      V            <- switch(uni.type,   constrained=P, single=1)
    }
    if(uni.prior   == "isotropic")   {
      psi.beta     <- matrix(vapply(Gseq, function(g) psi.beta[which.max(.ndeci(psi.beta[,g])),g], numeric(1L)), nrow=1, ncol=G)
    } else if(length(psi.beta) == 1) {
      psi.beta     <- matrix(psi.beta, nrow=1, ncol=G)
    }
    alpha.d1       <- cluster$alpha.d1
    alpha.d2       <- cluster$alpha.d2
    ad1.x          <- length(unique(alpha.d1)) == 1
    ad2.x          <- length(unique(alpha.d2)) == 1
    mu0g           <- cluster$l.switch[1]
    psi0g          <- cluster$l.switch[2]
    delta0g        <- cluster$l.switch[3]
    qstar0g        <- cluster$l.switch[4]
    label.switch   <- any(cluster$l.switch)
    eta            <- .sim_eta_p(N=N, Q=Q)
    phi            <- replicate(G, .sim_phi_p(Q=Q, P=P, nu=nu, plus1=nuplus1), simplify=FALSE)
    delta          <- lapply(Gseq, function(g) c(.sim_delta_p(alpha=alpha.d1[g], beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2[g], beta=beta.d2)))
    tau            <- lapply(delta, cumprod)
    lmat           <- lapply(Gseq, function(g) matrix(vapply(Pseq, function(j) .sim_load_ps(Q=Q, phi=phi[[g]][j,], tau=tau[[g]]), numeric(Q)), nrow=P, byrow=TRUE))
    psi.inv        <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    psi.inv        <- if(uni)     t(psi.inv)   else psi.inv
    psi.beta       <- if(one.uni) psi.beta[,1] else psi.beta
    if(Q < min(N - 1, Ledermann(P)))     {
      fact.ind     <- nn    <= P
      for(g in which(!fact.ind)) {
        fact       <- try(stats::factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]      <- fact$scores
          lmat[[g]]         <- unclass(fact$loadings)
          psi.inv[,g]       <- 1/fact$uniquenesses
        }
      }
    } else     {
      psi.tmp      <- psi.inv
      if(isTRUE(one.uni)) {
        psi.inv[,] <- 1/switch(uni.type, constrained=.col_vars(data), exp(mean(log(.col_vars(data)))))
      } else   {
        tmp.psi    <- ((nn[nn0] - 1)/(rowsum(data^2,  z) - rowsum(data, z)^2/nn[nn0]))
        psi.inv[,nn > 1]    <- tmp.psi[!is.nan(tmp.psi)]
      }
      inf.ind      <- is.infinite(psi.inv) | is.nan(psi.inv)
      psi.inv[inf.ind]      <- psi.tmp[inf.ind]
    }
    psi.inv[psi.inv == 0]   <- colMaxs(psi.inv[,which(psi.inv == 0, arr.ind=TRUE)[,2], drop=FALSE], value=TRUE)
    if(burnin       < 1)  {
      if(sw["mu.sw"])  mu.store[,,1]    <- mu
      if(sw["s.sw"])   eta.store[,,1]   <- eta
      if(sw["l.sw"])   load.store[,,,1] <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
      if(sw["psi.sw"]) psi.store[,,1]   <- 1/psi.inv
      if(sw["pi.sw"])  pi.store[,1]     <- pi.prop
      z.store[1,]           <- z
      Q0                    <- Qs > 0
      sigma                 <- if(uni) lapply(Gseq, function(g) as.matrix(1/psi.inv[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0)) else lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g]))
      log.probs             <- if(uni) vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log(pi.prop[g]), numeric(N)) else vapply(Gseq, function(g, Q=Q0[g]) { sigma <- if(Q) sigma[[g]] else sq_mat(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q) + log(pi.prop[g]) }, numeric(N))
      ll.store[1]           <- sum(rowLogSumExps(log.probs))
      Q.store[,1]           <- Qs
    }
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose   && iter   < burnin) utils::setTxtProgressBar(pb, iter)
      storage      <- is.element(iter, iters)

    # Mixing Proportions
      pi.prop      <- if(equal.pro) pi.prop else rDirichlet(G=G, alpha=pi.alpha, nn=nn)

    # Cluster Labels
      psi          <- 1/psi.inv
      Q0           <- Qs  > 0
      Q1           <- Qs == 1
      sigma        <- if(uni) lapply(Gseq, function(g) as.matrix(psi[,g] + if(Q0[g]) tcrossprod(lmat[[g]]) else 0)) else lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      log.pis      <- if(equal.pro) log.pis else log(pi.prop)
      if(uni) {
        log.probs  <- vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log.pis[g], numeric(N))
      } else  {
        log.probs  <- try(vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sq_mat(sigma[[g]]), log=TRUE, isChol=!Q) + log.pis[g], numeric(N)), silent=TRUE)
      }
      if(inherits(log.probs, "try-error")) {
        log.probs  <- vapply(Gseq, function(g, Q=Q0[g]) { sigma <- if(Q) sigma[[g]] else sq_mat(sigma[[g]]); dmvn(data, mu[,g], is.posi_def(sigma, make=TRUE)$X.new, log=TRUE, isChol=!Q) + log.pis[g] }, numeric(N))
      }
      z            <- gumbel_max(probs=log.probs)
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn  > 0
      nn.ind       <- which(nn0)
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-", check.margin=FALSE))
      if(!any(Q0))    {
        eta        <- .empty_mat(nr=N)
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
        lmat       <- replicate(G, .empty_mat(nr=P))
      } else {
        eta.tmp    <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else matrix(0, nrow=ifelse(Q0[g], 0, nn[g]), ncol=Qs[g]))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- lapply(Gseq, function(g) matrix(if(all(nn0[g], Q0[g])) vapply(Pseq, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j],
                             Q1=Q1[g], EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]]), numeric(Qs[g])) else
                             vapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]]), numeric(Qs[g])), nrow=P, byrow=TRUE))
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat      <- lapply(Gseq, function(g) { S <- c.data[[g]] - if(Q0[g]) tcrossprod(eta.tmp[[g]], lmat[[g]]) else 0; S * S } )
        psi.inv[,] <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,] <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta[,g], lmat=lmat[[g]],
                      P=P, Q0=Q0[g], eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
      }

    # Means
      sum.data     <- vapply(dat.g, colSums2, numeric(P))
      sum.data     <- if(uni) t(sum.data) else sum.data
      sum.eta      <- lapply(eta.tmp, colSums2)
      mu[,]        <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.sigma=mu.sigma, psi.inv=psi.inv[,g], mu.zero=mu.zero[,g], sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                             sum.data=sum.data[,g], lmat=lmat[[g]], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero[,g]), numeric(P))

    # Shrinkage
      if(all(Q0))     {
        load.2     <- lapply(lmat, .power2)
        phi        <- lapply(Gseq, function(g) if(nn0[g]) .sim_phi(Q=Qs[g], P=P, nu=nu, plus1=nuplus1,
                      tau=tau[[g]], load.2=load.2[[g]]) else .sim_phi_p(Q=Qs[g], P=P, nu=nu, plus1=nuplus1))

        sum.terms  <- lapply(Gseq, function(g) colSums2(phi[[g]] * load.2[[g]]))
        for(g in Gseq)  {
          Qg       <- Qs[g]
          Q1g      <- Q1[g]
          if(nn0[g])    {
            for(k in seq_len(Qg)) {
              delta[[g]][k] <- if(k > 1) .sim_deltak(alpha.d2=alpha.d2[g], beta.d2=beta.d2, delta.k=delta[[g]][k], tau.kq=tau[[g]][k:Qg], P=P,
                               Q=Qg, k=k, sum.term.kq=sum.terms[[g]][k:Qg]) else .sim_delta1(Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]],
                               alpha.d1=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), beta.d1=ifelse(Q1g, beta.d2, beta.d1), delta.1=delta[[g]][1])
              tau[[g]]      <- cumprod(delta[[g]])
            }
          } else {
            for(k in seq_len(Qg)) {
              delta[[g]][k] <- if(k > 1) .sim_delta_p(alpha=alpha.d2[g], beta=beta.d2) else .sim_delta_p(alpha=ifelse(Q1g, alpha.d2[g], alpha.d1[g]), beta=ifelse(Q1g, beta.d2, beta.d1))
              tau[[g]]      <- cumprod(delta[[g]])
            }
          }
        }
      }

    # Adaptation
      if(all(adapt, iter   > adaptat)) {
        if(stats::runif(1) < ifelse(iter < burnin, 0.5, exp(-b0 - b1 * (iter - adaptat)))) {
          colvec   <- lapply(nn.ind, function(g) (if(Q0[g]) colSums(abs(lmat[[g]])   < epsilon)/P else stats::runif(1)) >= prop)
          nonred   <- lapply(colvec, .which0)
          numred   <- lengths(colvec)  - lengths(nonred)
          notred   <- numred == 0
          ng.ind   <- seq_along(nn.ind)
          Qs.old   <- Qs[nn0]
          Qs[nn0]  <- pmax.int(0, vapply(ng.ind, function(h) if(notred[h]) Qs.old[h] + 1 else Qs.old[h] - numred[h], numeric(1L)))
          Q.big    <- Qs[nn0] > Q.star
          if((Q.bigs        <-  any(Q.big))) {
            notred <- notred & !Q.big
            Qs[nn0][Q.big]  <- Q.star
          }
          phi[nn0]          <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(phi[[g]][,seq_len(Qs.old[h])],  stats::rgamma(n=P, shape=nu + nuplus1, rate=nu)) else phi[[g]][,nonred[[h]], drop=FALSE])
          delta[nn0]        <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) c(delta[[g]][seq_len(Qs.old[h])],     stats::rgamma(n=1, shape=alpha.d2, rate=beta.d2)) else delta[[g]][nonred[[h]]])
          tau[nn0]          <- lapply(delta[nn.ind], cumprod)
          lmat[nn0]         <- lapply(nn.ind, function(g, h=which(nn.ind == g)) if(notred[h]) cbind(lmat[[g]][,seq_len(Qs.old[h])], stats::rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[h]], drop=FALSE])
          Qemp     <- Qs[!nn0]
          Qpop     <- Qs[nn0]
          Qmax     <- ifelse(all(Q.big), max(Qpop), max(Qpop[!Q.big]))
          Qmaxold  <- max(Qs.old)
         store.eta <- all(sw["s.sw"], storage)
          if(any(!nn0)  && Qmax  != max(Qemp))  {
           Qmaxseq <- seq_len(Qmax)
            for(g  in Gseq[!nn0][Qemp != Qmax]) {
              Qg   <- Qs[g]
              if(Qg > Qmax)  {
                phi[[g]]    <- phi[[g]][,Qmaxseq,  drop=FALSE]
                delta[[g]]  <- delta[[g]][Qmaxseq]
                tau[[g]]    <- tau[[g]][Qmaxseq]
                lmat[[g]]   <- lmat[[g]][,Qmaxseq, drop=FALSE]
              } else {
                while(Qg    != Qmax)   {
                 phi[[g]]   <- cbind(phi[[g]],  stats::rgamma(n=P, shape=nu + nuplus1,  rate=nu))
                 delta[[g]] <- c(delta[[g]],    stats::rgamma(n=1, shape=alpha.d2, rate=beta.d2))
                 tau[[g]]   <- cumprod(delta[[g]])
                 if(store.eta)         {
                  eta.tmp[[g]]   <- cbind(eta.tmp[[g]], .empty_mat(nc=1))
                 }
                 Qg         <- Qg + 1
                 lmat[[g]]  <- cbind(lmat[[g]], stats::rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qg] * tau[[g]][Qg]))))
                }
              }
            }
           Qs[Qmax  != Qs    & !nn0]  <- Qmax
          }
          if(store.eta)  {
            eta.tmp <- lapply(Gseq,  function(g) if(nn0[g] && Qs[g] > Qs.old[nn.ind == g]) cbind(eta.tmp[[g]], stats::rnorm(nn[g])) else eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE])
          }
        }
      }

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
          mu[,left]         <- mu[,right, drop=FALSE]
          lmat[left]        <- lmat[right]
          delta[left]       <- delta[right]
          phi[left]         <- phi[right]
          tau[left]         <- tau[right]
          psi.inv[,left]    <- psi.inv[,right, drop=FALSE]
          pi.prop[left]     <- pi.prop[right]
          nn[left]          <- nn[right]
          if(all(sw["s.sw"], storage)) {
            eta.tmp[left]   <- eta.tmp[right]
            dat.g[left]     <- dat.g[right]
            Q0[left]        <- Q0[right]
            nn0[left]       <- nn0[right]
          }
          if(mu0g)         {
           mu.zero[,left]   <- mu.zero[,right, drop=FALSE]
          }
          if(psi0g)        {
           psi.beta[,left]  <- psi.beta[,right, drop=FALSE]
          }
          if(all(delta0g,
                 !ad1.x))  {
           alpha.d1[left]   <- alpha.d1[right]
          }
          if(all(delta0g,
                 !ad2.x))  {
           alpha.d2[left]   <- alpha.d2[right]
          }
        }
      }

      if(Q.bigs && !Q.large && iter > burnin) {        warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large    <- TRUE
      }
      if(z.err  && !err.z) {                           warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z      <- TRUE
      }
      if(storage)  {
        if(verbose)   utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])        mu.store[,,new.it]   <- mu
        if(all(sw["s.sw"],
           any(Q0)))     {
          max.Q    <- max(Qs)
          eta.tmp  <- if(length(unique(Qs))  != 1)     lapply(Gseq,       function(g) cbind(eta.tmp[[g]], matrix(0, nrow=nn[g], ncol=max.Q - Qs[g]))) else eta.tmp
          q0ng     <- (!Q0   | Qs[Gseq] == 0) & nn0[Gseq]
          if(any(q0ng))  {
            eta.tmp[q0ng]   <-                         lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                eta.store[,seq_len(max.Q),new.it]   <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])   {
          for(g in Gseq) {
            Qseqg  <- seq_len(Qs[g])
            if(Q0[g]) load.store[,Qseqg,g,new.it]   <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])   psi.store[,,new.it]      <- 1/psi.inv
        if(sw["pi.sw"])    pi.store[,new.it]        <- pi.prop
                           z.store[new.it,]         <- as.integer(z)
                           ll.store[new.it]         <- sum(rowLogSumExps(log.probs))
                           Q.store[,new.it]         <- as.integer(Qs)
      }
    }
    if(verbose)       close(pb)
    Qmax           <- seq_len(max(Q.store))
    eta.store      <- if(sw["s.sw"])  tryCatch(eta.store[,Qmax,, drop=FALSE],   error=function(e) eta.store)
    load.store     <- if(sw["l.sw"])  tryCatch(load.store[,Qmax,,, drop=FALSE], error=function(e) load.store)
    returns        <- list(mu       = if(sw["mu.sw"])  tryCatch(provideDimnames(mu.store,  base=list(varnames, "", ""), unique=FALSE), error=function(e) mu.store),
                           eta      = if(sw["s.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(eta.store),            error=function(e) eta.store),  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                           load     = if(sw["l.sw"])   tryCatch(provideDimnames(tryCatch(as.simple_sparse_array(load.store),           error=function(e) load.store), base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                           psi      = if(sw["psi.sw"]) tryCatch(provideDimnames(psi.store, base=list(varnames, "", ""), unique=FALSE), error=function(e) psi.store),
                           pi.prop  = if(sw["pi.sw"])  pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           Q.store  = Q.store,
                           time     = init.time)
    attr(returns, "Q.big")  <- Q.large
    return(returns)
  }
