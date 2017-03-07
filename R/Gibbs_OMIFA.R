###########################################################################
### Gibbs Sampler for Overfitted Mixtures of Infinite Factor Analysers ####
###########################################################################

# Gibbs Sampler Function
  .gibbs_OMIFA       <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu, uni.type,
                                 uni.prior, burnin, thinning, adapt, psi.alpha, psi.beta,
                                 verbose, alpha.d1, alpha.d2, sw, cluster, nu, b0, b1, mu,
                                 prop, beta.d1, beta.d2, adapt.at, epsilon, nuplus1, ...) {

  # Define & initialise variables
    start.time       <- proc.time()
    total            <- max(iters)
    if(verbose)         pb    <- txtProgressBar(min=0, max=total, style=3)
    n.store          <- length(iters)
    Gseq             <- seq_len(G)
    Pseq             <- seq_len(P)
    Nseq             <- seq_len(N)
    obsnames         <- rownames(data)
    varnames         <- colnames(data)
    colnames(data)   <- NULL
    if(sw["mu.sw"])  {
      mu.store       <- array(0, dim=c(P, G, n.store))
    }
    if(sw["s.sw"])   {
      eta.store      <- array(0, dim=c(N, Q, n.store))
    }
    if(sw["l.sw"])   {
      load.store     <- array(0, dim=c(P, Q, G, n.store))
    }
    if(sw["psi.sw"]) {
      psi.store      <- array(0, dim=c(P, G, n.store))
    }
    if(sw["pi.sw"])  {
      pi.store       <- matrix(0, nrow=G, ncol=n.store)
    }
    z.store          <- matrix(0, nrow=N, ncol=n.store)
    ll.store         <- rep(0, n.store)
    Q.star           <- Q
    Qs               <- rep(Q, G)
    Q.store          <- matrix(0, nrow=G, ncol=n.store)
    Q.large          <- Q.big <- Q.bigs <- FALSE
    err.z            <- z.err <- FALSE
    G.store          <- ll.store

    mu.sigma         <- 1/sigma.mu
    sig.mu.sqrt      <- sqrt(sigma.mu)
    z                <- cluster$z
    z.temp           <- factor(z, levels=Gseq)
    nn               <- tabulate(z, nbins=G)
    nn.ind           <- which(nn > 0)
    pi.prop          <- cluster$pi.prop
    pi.alpha         <- cluster$pi.alpha
    .sim_psi_inv     <- switch(uni.type,  unconstrained=.sim_psi_iu,  isotropic=.sim_psi_ii)
    .sim_psi_ip      <- switch(uni.type,  unconstrained=.sim_psi_ipu, isotropic=.sim_psi_ipi)
    psi.beta         <- switch(uni.prior, isotropic=unique(round(psi.beta, min(nchar(psi.beta)))), psi.beta)
    eta              <- .sim_eta_p(N=N, Q=Q)
    phi              <- lapply(Gseq, function(g) .sim_phi_p(Q=Q, P=P, nu=nu, plus1=nuplus1))
    delta            <- lapply(Gseq, function(g) c(.sim_delta_p(alpha=alpha.d1, beta=beta.d1), .sim_delta_p(Q=Q, alpha=alpha.d2, beta=beta.d2)))
    tau              <- lapply(delta, cumprod)
    lmat             <- lapply(Gseq, function(g) matrix(unlist(lapply(Pseq, function(j) .sim_load_ps(Q=Q, phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nrow=P, byrow=TRUE))
    psi.inv          <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q < P - sqrt(P + Q))   {
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
      psi.inv        <- vapply(Gseq, function(g) if(nn[g] > 1) 1/Rfast::colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
      inf.ind        <- is.infinite(psi.inv)
      psi.inv[inf.ind]        <- psi.tmp[inf.ind]
    }
    index          <- order(nn, decreasing=TRUE)
    pi.prop        <- pi.prop[index]
    mu             <- mu[,index, drop=FALSE]
    phi            <- phi[index]
    delta          <- delta[index]
    tau            <- tau[index]
    lmat           <- lmat[index]
    Qs             <- Qs[index]
    psi.inv        <- psi.inv[,index, drop=FALSE]
    z              <- factor(z, labels=match(nn.ind, index))
    z              <- as.numeric(levels(z))[z]
    G.non          <- G
    if(burnin       < 1)  {
      mu.store[,,1]           <- mu
      eta.store[,,1]          <- eta
      load.store[,,,1]        <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
      psi.store[,,1]          <- 1/psi.inv
      pi.store[,1]            <- pi.prop
      z.store[,1]             <- z
      sigma                   <- lapply(Gseq, function(g) make.positive.definite(tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])))
      Q0                      <- Qs > 0
      log.probs               <- vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
      ll.store[1]             <- sum(gumbel_max(probs=log.probs, log.like=TRUE)$log.like)
      Q.store[,1]             <- Qs
      G.store[1]              <- G.non
    }
    init.time        <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total)[-1]) {
      if(verbose     && iter   < burnin) setTxtProgressBar(pb, iter)
      storage        <- is.element(iter, iters)

    # Mixing Proportions & Re-ordering
      pi.prop        <- if(G  == 1) 1 else rDirichlet(G=G, alpha=pi.alpha, nn=nn)
      index          <- order(nn, decreasing=TRUE)
      pi.prop        <- pi.prop[index]
      mu             <- mu[,index, drop=FALSE]
      phi            <- phi[index]
      delta          <- delta[index]
      tau            <- tau[index]
      lmat           <- lmat[index]
      Qs             <- Qs[index]
      psi.inv        <- psi.inv[,index, drop=FALSE]
      z              <- factor(z, labels=match(nn.ind, index))
      z              <- as.numeric(levels(z))[z]

    # Cluster Labels
      psi            <- 1/psi.inv
      Q0             <- Qs  > 0
      Q1             <- Qs == 1
      if(G > 1)   {
        sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
        log.check    <- capture.output(log.probs <- try(vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) sigma[[g]] else sqrt(sigma[[g]]), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N)), silent=TRUE))
        if(inherits(log.probs, "try-error")) {
          log.probs  <- vapply(Gseq, function(g, Q=Q0[g]) dmvn(data, mu[,g], if(Q) make.positive.definite(sigma[[g]]) else make.positive.definite(sqrt(sigma[[g]])), log=TRUE, isChol=!Q) + log(pi.prop[g]), numeric(N))
        }
        z.res        <- gumbel_max(probs=log.probs, log.like=TRUE)
        z            <- z.res$z
      } else      {
        z            <- rep(1, N)
      }
      nn             <- tabulate(z, nbins=G)
      nn0            <- nn  > 0
      nn.ind         <- which(nn0)
      G.non          <- length(nn.ind)
      dat.g          <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])

    # Scores & Loadings
      c.data         <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta          <- base::matrix(0, nrow=N, ncol=0)
        eta.tmp      <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
        lmat         <- lapply(Gseq, base::matrix, 0, nrow=P, ncol=0)
      } else {
        eta.tmp      <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) .sim_score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else base::matrix(0, nrow=ifelse(Q0[g], 0, nn[g]), ncol=Qs[g]))
        EtE          <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat         <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Pseq, function(j) .sim_load_s(Q=Qs[g], c.data=c.data[[g]][,j], Q1=Q1[g],
                               EtE=EtE[[g]], eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nrow=P, byrow=TRUE) else
                               base::matrix(unlist(lapply(Pseq, function(j) .sim_load_ps(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]])), use.names=FALSE), nrow=P, byrow=FALSE))
      }

    # Uniquenesses
      psi.inv        <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                               P=P, eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))

    # Means
      sum.data       <- vapply(dat.g, colSums, numeric(P))
      sum.eta        <- lapply(eta.tmp, colSums)
      mu             <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                               sum.data=sum.data[,g], lmat=lmat[[g]], mu.zero=mu.zero) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero), numeric(P))

    # Local Shrinkage
      load.2         <- lapply(lmat, .power2)
      phi            <- lapply(Gseq, function(g) if(nn0[g]) .sim_phi(Q=Qs[g], P=P, nu=nu, plus1=nuplus1,
                        tau=tau[[g]], load.2=load.2[[g]]) else .sim_phi_p(Q=Qs[g], P=P, nu=nu, plus1=nuplus1))

    # Global Shrinkage
      sum.terms      <- lapply(Gseq, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gseq) {
        Qg           <- Qs[g]
        Q1g          <- Q1[g]
        if(nn0[g])   {
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
            Qs[Qmax   < Qs & !nn0]  <- Qmax
            Qmaxseq  <- seq_len(Qmax)
            for(g in Gseq[!nn0][Qemp > Qmax]) {
              phi[[g]]        <- phi[[g]][,Qmaxseq,  drop=FALSE]
              delta[[g]]      <- delta[[g]][Qmaxseq]
              tau[[g]]        <- tau[[g]][Qmaxseq]
              lmat[[g]]       <- lmat[[g]][,Qmaxseq, drop=FALSE]
            }
          }
          if(all(sw["s.sw"], storage)) {
            eta.tmp  <- lapply(Gseq,   function(g) if(nn0[g] && Qs[g] > Qs.old[which(nn.ind == g)]) cbind(eta.tmp[[g]], rnorm(nn[g])) else eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE])
          }
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
           any(Q0)))     {
          max.Q      <-  max(Qs)
          eta.tmp    <-  if(length(unique(Qs)) != 1)      lapply(Gseq,       function(g) cbind(eta.tmp[[g]], base::matrix(0, nrow=nn[g], ncol=max.Q - Qs[g]))) else eta.tmp
          q0ng       <-  (!Q0  | Qs[Gseq] == 0) & nn0[Gseq]
          if(any(q0ng))  {
            eta.tmp[q0ng]     <-                          lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- row.names(dat.g[[g]]); x })
          }
                   eta.store[,seq_len(max.Q),new.it]   <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
        }
        if(sw["l.sw"])   {
          for(g in Gseq) {
            Qseqg    <-  seq_len(Qs[g])
            if(Q0[g])    load.store[,Qseqg,g,new.it]   <- lmat[[g]]
          }
        }
        if(sw["psi.sw"]) psi.store[,,new.it]           <- 1/psi.inv
        if(sw["pi.sw"])  pi.store[,new.it]             <- pi.prop
                         z.store[,new.it]              <- z
                         ll.store[new.it]              <- sum(z.res$log.like)
                         Q.store[,new.it]              <- Qs
                         G.store[new.it]               <- G.non
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
                             z.store   = z.store,
                             ll.store  = ll.store,
                             G.store   = G.store,
                             Q.store   = tryCatch(Q.store[Gmax,, drop=FALSE],          error=function(e) Q.store),
                             time      = init.time)
    attr(returns, "Q.big")    <- Q.large
    return(returns)
  }
