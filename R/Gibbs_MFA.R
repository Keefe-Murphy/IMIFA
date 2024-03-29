#########################################################################
### Gibbs Sampler for (Finite) Mixtures of (Finite) Factor Analysers ####
#########################################################################

# Gibbs Sampler Function
  .gibbs_MFA       <- function(Q, data, iters, N, P, G, mu.zero, sigma.mu, mu,
                               sigma.l, burnin, thinning, uni.type, uni.prior, col.mean,
                               sw, psi.alpha, psi.beta, verbose, cluster, equal.pro, ...) {

  # Define & initialise variables
    start.time     <- proc.time()
    sq_mat         <- if(P  > 50) function(x) diag(sqrt(diag(x))) else sqrt
    matrix         <- base::matrix
    total          <- max(iters)
    if(verbose)       pb   <- utils::txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    Nseq           <- seq_len(N)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    Q0             <- Q  > 0
    Q1             <- Q == 1
    uni            <- P == 1
    sw["s.sw"]     <- sw["s.sw"] && Q0
    sw["l.sw"]     <- sw["l.sw"] && Q0
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
    err.z          <- zerr <- FALSE
    ll.store       <- integer(n.store)

    mu.sigma       <- 1/sigma.mu
    sig.mu.sqrt    <- sqrt(sigma.mu)
    l.sigma        <- diag(1/sigma.l, Q)
    sig.l.sqrt     <- sqrt(sigma.l)
    if(all(mu.zero  == 0))  {
      mu.zero      <- matrix(0L, nrow=1L, ncol=G)
      cluster$l.switch[1L] <- FALSE
    }
    if(length(mu.zero)  == 1) {
      mu.zero      <- matrix(mu.zero, nrow=1L, ncol=G)
    }
    mu.prior       <- mu.sigma * mu.zero
    z              <- cluster$z
    nn             <- tabulate(z, nbins=G)
    nn0            <- nn > 0
    z.temp         <- factor(z, levels=Gseq)
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
    mu0g           <- cluster$l.switch[1L]
    psi0g          <- cluster$l.switch[2L]
    label.switch   <- any(cluster$l.switch)
    if(Q0)                  {
      eta          <- .sim_eta_p(N=N, Q=Q)
      lmat         <- array(vapply(Gseq, function(g) .sim_load_p(Q=Q, P=P, sig.l.sqrt=sig.l.sqrt), numeric(P * Q)), dim=c(P, Q, G))
    } else                  {
      eta          <- .empty_mat(nr=N)
      lmat         <- array(0L, dim=c(P, 0, G))
    }
    if(isTRUE(one.uni))     {
      psi.beta     <- psi.beta[,1L]
      psi.inv      <- matrix(.sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), nrow=P, ncol=G)
    } else psi.inv <- vapply(Gseq, function(g) .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    psi.inv        <- if(uni)     t(psi.inv)    else psi.inv
    if(isTRUE(one.uni))     {
      psi.inv[]    <- 1/switch(EXPR=uni.type, constrained=colVars(data, center=col.mean, refine=FALSE, useNames=FALSE), max(colVars(data, center=col.mean, refine=FALSE, useNames=FALSE)))
    } else  {
      tmp.psi      <- (nn[nn0] - 1L)/pmax(rowsum(data^2, z) - rowsum(data, z)^2/nn[nn0], 0L)
      tmp.psi      <- switch(EXPR=uni.type, unconstrained=t(tmp.psi), matrix(Rfast::rowMaxs(tmp.psi, value=TRUE), nrow=P, ncol=G, byrow=TRUE))
      psi.inv[,nn   > 1]   <- tmp.psi[!is.nan(tmp.psi)]
      rm(tmp.psi)
    }
    max.p          <- (psi.alpha  - 1)/psi.beta
    inf.ind        <- psi.inv > max(max.p)
    psi.inv[inf.ind]       <- matrix(max.p, nrow=P, ncol=G)[inf.ind]
    rm(max.p, inf.ind)
    init.time      <- proc.time() - start.time

  # Iterate
    for(iter in seq_len(total))   {
      if(verbose   && iter  < burnin)  utils::setTxtProgressBar(pb, iter)
      storage      <- is.element(iter, iters)

    # Mixing Proportions
      pi.prop      <- if(equal.pro) pi.prop else rDirichlet(G=G, alpha=pi.alpha, nn=nn)

    # Scores & Loadings
      dat.g        <- lapply(Gseq, function(g) data[z == g,, drop=FALSE])
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2L, mu[,g], FUN="-", check.margin=FALSE))
      if(Q0) {
        eta.tmp    <- lapply(Gseq, function(g) if(nn0[g]) .sim_score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else .empty_mat(nc=Q))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) matrix(if(nn0[g]) vapply(Pseq, function(j) .sim_load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]], Q1=Q1,
                      EtE=EtE[[g]], psi.inv=psi.inv[,g][j]), numeric(Q)) else .sim_load_p(Q=Q, P=P, sig.l.sqrt=sig.l.sqrt), nrow=P, byrow=TRUE)), use.names=FALSE), dim=c(P, Q, G))
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else {
        eta.tmp    <- lapply(Gseq, function(g) eta[z == g,, drop=FALSE])
      }

    # Uniquenesses
      if(isTRUE(one.uni)) {
        S.mat      <- lapply(Gseq, function(g) { S   <- c.data[[g]] - if(Q0) tcrossprod(eta.tmp[[g]], if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else 0L; S^2 } )
        psi.inv[,] <- .sim_psi_inv(uni.shape, psi.beta, S.mat, V)
      } else {
        psi.inv[,] <- vapply(Gseq, function(g) if(nn0[g]) .sim_psi_inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], eta=eta.tmp[[g]], psi.beta=psi.beta[,g],
                      P=P, Q0=Q0, lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim_psi_ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
      }

    # Means
      sum.data     <- vapply(dat.g, colSums2, useNames=FALSE, numeric(P))
      sum.data     <- if(uni) t(sum.data) else sum.data
      sum.eta      <- lapply(eta.tmp, colSums2, useNames=FALSE)
      mu[,]        <- vapply(Gseq, function(g) if(nn0[g]) .sim_mu(mu.prior=mu.prior[,g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                             lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], N=nn[g], P=P) else .sim_mu_p(P=P, sig.mu.sqrt=sig.mu.sqrt, mu.zero=mu.zero[,g]), numeric(P))

    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- if(uni) lapply(Gseq, function(g) as.matrix(psi[,g] + if(Q0) tcrossprod(as.matrix(lmat[,,g])) else 0L)) else lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      log.pis      <- if(equal.pro) log.pis else log(pi.prop)
      if(uni) {
        log.probs  <- vapply(Gseq, function(g) stats::dnorm(data, mu[,g], sq_mat(sigma[[g]]), log=TRUE) + log.pis[g], numeric(N))
      } else  {
        log.probs  <- try(vapply(Gseq, function(g) dmvn(data, mu[,g], if(Q0) sigma[[g]] else sq_mat(sigma[[g]]), log=TRUE, isChol=!Q0) + log.pis[g], numeric(N)), silent=TRUE)
        if(zerr    <- inherits(log.probs, "try-error")) {
         log.probs <- vapply(Gseq, function(g) { sigma <- if(Q0) is.posi_def(sigma[[g]], make=TRUE)$X.new else sq_mat(sigma[[g]]); dmvn(data, mu[,g], sigma, log=TRUE, isChol=!Q0) + log.pis[g] }, numeric(N))
        }
      }
      z            <- gumbel_max(probs=log.probs)

    # Label Switching
      if(label.switch) {
        sw.lab     <- .lab_switch(z.new=z, z.old=z.temp)
        z.perm     <- sw.lab$z.perm
        left       <- as.integer(unname(z.perm))
        right      <- as.integer(names(z.perm))
        if(!identical(left, right)) {
          z        <- sw.lab$z
          mu[,left]        <- mu[,right,       drop=FALSE]
          lmat[,,left]     <- lmat[,,right,    drop=FALSE]
          psi.inv[,left]   <- psi.inv[,right,  drop=FALSE]
          pi.prop[left]    <- pi.prop[right]
          if(mu0g)  {
           mu.zero[,left]  <- mu.zero[,right,  drop=FALSE]
           mu.prior[,left] <- mu.prior[,right, drop=FALSE]
          }
          if(psi0g) {
           psi.beta[,left] <- psi.beta[,right, drop=FALSE]
          }
        }
      }
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0

      if(zerr && !err.z) {                                    cat("\n"); warning("\nAlgorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices\n", call.=FALSE, immediate.=TRUE)
        err.z      <- TRUE
      }
      if(storage)  {
        if(verbose)   utils::setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])               mu.store[,,new.it]   <- mu
        if(all(sw["s.sw"], Q0))      eta.store[,,new.it]   <- eta
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])             psi.store[,,new.it]   <- 1/psi.inv
        if(sw["pi.sw"])                pi.store[,new.it]   <- pi.prop
                                        z.store[new.it,]   <- as.integer(z)
                                        ll.store[new.it]   <- sum(rowLogSumExps(log.probs, useNames=FALSE))
      }
    }
    if(verbose)       close(pb)
    returns        <- list(mu       = if(sw["mu.sw"])         tryCatch(provideDimnames(mu.store,   base=list(varnames, "", ""),     unique=FALSE), error=function(e) mu.store),
                           eta      = if(all(sw["s.sw"], Q0)) tryCatch(provideDimnames(eta.store,  base=list(obsnames, "", ""),     unique=FALSE), error=function(e) eta.store),
                           load     = if(all(sw["l.sw"], Q0)) tryCatch(provideDimnames(load.store, base=list(varnames, "", "", ""), unique=FALSE), error=function(e) load.store),
                           psi      = if(sw["psi.sw"])        tryCatch(provideDimnames(psi.store,  base=list(varnames, "", ""),     unique=FALSE), error=function(e) psi.store),
                           pi.prop  = if(sw["pi.sw"])         pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           time     = init.time)
    attr(returns, "K")  <- PGMM_dfree(Q=Q, P=P, G=G, method=switch(EXPR=uni.type, unconstrained="UUU", isotropic="UUC", constrained="UCU", single="UCC"), equal.pro=equal.pro)
      return(returns)
  }
