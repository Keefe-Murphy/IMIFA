#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  .gibbs_IMFA       <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, alpha.step, discount,
                                a.hyper, mu, sigma.mu, burnin, thinning, trunc.G, d.hyper, learn.d,
                                ind.slice, psi.alpha, psi.beta, verbose, sw, cluster, DP.lab.sw, ...) {
        
  # Define & initialise variables
    start.time      <- proc.time()
    total           <- max(iters)
    if(verbose)        pb   <- txtProgressBar(min=0, max=total, style=3)
    n.store         <- length(iters)
    Gs              <- seq_len(G)
    Ts              <- seq_len(trunc.G)
    Ps              <- seq_len(P)
    Ns              <- seq_len(N)
    obsnames        <- rownames(data)
    varnames        <- colnames(data)
    colnames(data)  <- NULL
    facnames        <- paste0("Factor ", seq_len(Q))
    gnames          <- paste0("Group ", Ts)
    iternames       <- paste0("Iteration", seq_len(n.store))
    Q0              <- Q  > 0
    Q0s             <- rep(Q0, trunc.G)
    Q1              <- Q == 1
    if(sw["mu.sw"])  {
      mu.store      <- provideDimnames(array(0, dim=c(P, trunc.G, n.store)), base=list(varnames, gnames, iternames))
    }
    if(sw["s.sw"])   {
      eta.store     <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)   <- list(obsnames, if(Q0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store    <- array(0, dim=c(P, Q, trunc.G, n.store))
      dimnames(load.store)  <- list(varnames, if(Q0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store     <- provideDimnames(array(0, dim=c(P, trunc.G, n.store)), base=list(varnames, gnames, iternames))
    }
    if(sw["pi.sw"])  {
      pi.store      <- provideDimnames(matrix(0, nr=trunc.G, nc=n.store), base=list(gnames, iternames))
    }
    z.store         <- provideDimnames(matrix(0, nr=N, nc=n.store), base=list(obsnames, iternames))
    ll.store        <- setNames(rep(0, n.store), iternames)
    acc1            <- acc2 <- FALSE
    err.z           <- zerr <- FALSE
    G.store         <- setNames(rep(0, n.store), iternames)
    not.fixed       <- alpha.step != "fixed"
    if(not.fixed) {
      alpha.store   <- setNames(rep(0, n.store), iternames)
      alpha.shape   <- a.hyper[1]
      alpha.rate    <- a.hyper[2]
    }
    if(learn.d)   {
      d.store       <- setNames(rep(0, n.store), iternames)
      d.shape1      <- d.hyper[1]
      d.shape2      <- d.hyper[2]
    }
    MH.step         <- alpha.step == "metropolis"
    if(MH.step)   {
      rate          <- setNames(rep(0, n.store), iternames)
    }
    if(DP.lab.sw) {
      lab.rate      <- matrix(0, nr=2, nc=n.store)  
    }
    mu.sigma        <- 1/sigma.mu
    z               <- cluster$z
    pi.alpha        <- cluster$pi.alpha
    pi.prop         <- c(cluster$pi.prop, .sim.pi(pi.alpha=pi.alpha, nn=rep(0, trunc.G), N=N, inf.G=TRUE, len=trunc.G, discount=discount)$pi.prop[-Gs])
    nn              <- tabulate(z, nbins=trunc.G)
    mu              <- cbind(mu, vapply(seq_len(trunc.G - G), function(g) .sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P)))
    eta             <- .sim.eta.p(N=N, Q=Q)
    lmat            <- lapply(Ts, function(t) .sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE))
    psi.inv         <- vapply(Ts, function(t) .sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    if(Q0 && Q   < P - sqrt(P + Q)) {
      for(g in which(nn      > P))  {
        fact        <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
        if(!inherits(fact, "try-error")) {
          eta[z == g,]      <- fact$scores
          lmat[[g]]         <- fact$loadings
          psi.inv[,g]       <- 1/fact$uniquenesses
        } 
      }
    } else {
      psi.tmp       <- psi.inv
      psi.inv[,Gs]  <- vapply(Gs, function(g) if(nn[g] > 1) 1/colVars(data[z == g,, drop=FALSE]) else psi.tmp[,g], numeric(P))
      inf.ind       <- is.infinite(psi.inv)
      psi.inv[inf.ind]      <- psi.tmp[inf.ind]
    }
    l.sigma         <- diag(1/sigma.l, Q)
    lmat            <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
    index           <- order(pi.prop, decreasing=TRUE)
    pi.prop         <- pi.prop[index]
    mu              <- mu[,index, drop=FALSE]
    lmat            <- lmat[,,index, drop=FALSE]
    psi.inv         <- psi.inv[,index, drop=FALSE]
    nn              <- nn[index]
    nn0             <- nn > 0
    nn.ind          <- which(nn0)
    G.non           <- sum(nn0)
    z               <- factor(z, labels=match(nn.ind, index))
    z               <- as.numeric(levels(z))[z]
    ksi             <- (1 - rho) * rho^(Ts - 1)
    log.ksi         <- log(ksi)
    slice.logs      <- c(- Inf, 0L)
    if(burnin        < 1)  {
      mu.store[,,1]         <- mu
      eta.store[,,1]        <- eta
      load.store[,,,1]      <- lmat
      psi.store[,,1]        <- 1/psi.inv
      pi.store[,1]          <- pi.prop/sum(pi.prop)
      z.store[,1]           <- z
      ll.store[1]           <- sum(.sim.z(data=data, mu=mu[,Gs], Gseq=Gs, N=N, pi.prop=pi.store[Gs,1], sigma=lapply(Gs, function(g)
                                         make.positive.definite(tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g]))), Q0=Q0s)$log.likes)
      G.store[1]            <- G.non
      if(not.fixed) {
        alpha.store[1]      <- pi.alpha 
      }
      if(learn.d)   {
        d.store[1]          <- discount
      }
    }
    init.time       <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose    && iter  < burnin) setTxtProgressBar(pb, iter)
    
    # Mixing Proportions
      weights       <- .sim.pi(pi.alpha=pi.alpha, nn=nn, N=N, inf.G=TRUE, len=trunc.G, discount=discount)
      pi.prop       <- weights$pi.prop
      Vs            <- weights$Vs
      
    # Re-ordering & Slice Sampler
      index         <- order(pi.prop, decreasing=TRUE)
      pi.prop       <- pi.prop[index]
      Vs            <- Vs[index]
      mu            <- mu[,index, drop=FALSE]
      lmat          <- lmat[,,index, drop=FALSE]
      psi.inv       <- psi.inv[,index, drop=FALSE]
      z             <- factor(z, labels=match(nn.ind, index))
      z             <- as.numeric(levels(z))[z]
      if(!ind.slice) {
        ksi         <- pi.prop
        log.ksi     <- log(ksi)
      }
      u.slice       <- runif(N, 0, ksi[z])
      G             <- max(vapply(Ns, function(i) sum(u.slice[i] < ksi), integer(1L)))
      Gs            <- seq_len(G)
      log.slice.ind <- vapply(Gs, function(g) slice.logs[1 + (u.slice < ksi[g])] - log.ksi[g], numeric(N))
    
    # Cluster Labels
      psi           <- 1/psi.inv
      sigma         <- lapply(Gs, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.log         <- capture.output({ z.res <- try(.sim.z(data=data, mu=mu[,Gs, drop=FALSE], sigma=sigma, Gseq=Gs, N=N, pi.prop=pi.prop[Gs], log.slice.ind=log.slice.ind, Q0=Q0s[Gs]), silent=TRUE) })
      zerr          <- inherits(z.res, "try-error")
      if(zerr) {
        sigma       <- lapply(sigma, make.positive.definite)
        z.res       <- .sim.z(data=data, mu=mu[,Gs, drop=FALSE], sigma=sigma, Gseq=Gs, N=N, pi.prop=pi.prop[Gs], log.slice.ind=log.slice.ind, Q0=Q0s[Gs])
      }
      z             <- z.res$z
      nn            <- tabulate(z, nbins=trunc.G)
      nn0           <- nn > 0
      nn.ind        <- which(nn0)
      dat.g         <- lapply(Gs, function(g) data[z == g,, drop=FALSE])
      
    # Scores & Loadings
      c.data        <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0)   {
        eta.tmp     <- lapply(Gs, function(g) if(nn0[g]) .sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else base::matrix(0, nr=0, nc=Q))
        EtE         <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[,,Gs]  <- array(unlist(lapply(Gs, function(g) if(nn0[g]) matrix(unlist(lapply(Ps, function(j) .sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], eta=eta.tmp[[g]], Q1=Q1, EtE=EtE[[g]], 
                             P=P, psi.inv=psi.inv[,g][j], shrink=FALSE)), use.names=FALSE), nr=P, byrow=TRUE) else .sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=FALSE)), use.names=FALSE), dim=c(P, Q, G))
        eta         <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      } else   {
        eta.tmp     <- lapply(Gs, function(g) eta[z == g,, drop=FALSE])
      }
      
    # Means
      sum.data      <- vapply(dat.g, colSums, numeric(P))
      sum.eta       <- lapply(eta.tmp, colSums)
      mu[,Gs]       <- vapply(Gs, function(g) if(nn0[g]) .sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[,g], sum.eta=sum.eta[[g]],
                              lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g], mu.zero=mu.zero) else .sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P))
                    
    # Uniquenesses
      psi.inv[,Gs]  <- vapply(Gs, function(g) if(nn0[g]) .sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, eta=eta.tmp[[g]], psi.beta=psi.beta,
                              lmat=if(Q1) as.matrix(lmat[,,g]) else lmat[,,g]) else .sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
      
    # Alpha
      if(not.fixed) {
        if(MH.step) {
          MH.alpha  <- .sim.alpha.m(alpha=pi.alpha, lower=alpha.shape, upper=alpha.rate, trunc.G=trunc.G, Vs=Vs, discount=discount) 
          pi.alpha  <- MH.alpha$alpha  
        } else {
          pi.alpha  <- .sim.alpha.g(alpha=pi.alpha, shape=alpha.shape, rate=alpha.rate, G=G, N=N, discount=discount) 
        }
      }
      
    # Label Switching
      G.non         <- sum(nn0)
      if(DP.lab.sw  && G.non > 1)  {
        move1       <- .lab.move1(nn.ind=nn.ind, pi.prop=pi.prop, nn=nn)
        acc1        <- move1$rate1
        if(acc1)    {
          sw1       <- move1$sw
          sw1x      <- c(sw1[2], sw1[1])
          nn[sw1]   <- nn[sw1x]
          nn0[sw1]  <- nn0[sw1x]
          nn.ind    <- which(nn0)
          Vs[sw1]   <- Vs[sw1x]
          mu[,sw1]  <- mu[,sw1x, drop=FALSE]
          lmat[,,sw1]       <- lmat[,,sw1x,   drop=FALSE]
          psi.inv[,sw1]     <- psi.inv[,sw1x, drop=FALSE]
          pi.prop[sw1]      <- pi.prop[sw1x]
          zsw1      <- z == sw1[1]
          z[z == sw1[2]]    <- sw1[1]
          z[zsw1]   <- sw1[2]            
        } 
        move2       <- .lab.move2(nn.ind=nn.ind, Vs=Vs, nn=nn)
        acc2        <- move2$rate2
        if(acc2)    {
          sw2       <- move2$sw
          sw2x      <- c(sw2[2], sw2[1])
          nn[sw2]   <- nn[sw2x]
          nn0[sw2]  <- nn0[sw2x]
          nn.ind    <- which(nn0)
          mu[,sw2]  <- mu[,sw2x, drop=FALSE]
          lmat[,,sw2]       <- lmat[,,sw2x,   drop=FALSE]
          psi.inv[,sw2]     <- psi.inv[,sw2x, drop=FALSE]
          pi.prop[sw2]      <- pi.prop[sw2x]
          zsw2      <- z == sw2[1]
          z[z == sw2[2]]    <- sw2[1]
          z[zsw2]   <- sw2[2]            
        }
      }
      if(zerr && !err.z) {     warning("Algorithm may slow due to corrections for Choleski decompositions of non-positive-definite covariance matrices", call.=FALSE)
        err.z       <- TRUE
      }
      
      if(is.element(iter, iters)) {
        if(verbose)    setTxtProgressBar(pb, iter)
        new.it      <- which(iters == iter)
        if(sw["mu.sw"])               mu.store[,,new.it]    <- mu 
        if(all(sw["s.sw"], Q0))       eta.store[,,new.it]   <- eta
        if(all(sw["l.sw"], Q0))       load.store[,,,new.it] <- lmat
        if(sw["psi.sw"])              psi.store[,,new.it]   <- psi
        if(sw["pi.sw"])               pi.store[,new.it]     <- pi.prop
        if(not.fixed)                 alpha.store[new.it]   <- pi.alpha
        if(learn.d)                   d.store[new.it]       <- discount
        if(MH.step)                   rate[new.it]          <- MH.alpha$rate
        if(DP.lab.sw)                 lab.rate[,new.it]     <- c(acc1, acc2)
                                      z.store[,new.it]      <- z 
                                      ll.store[new.it]      <- sum(z.res$log.likes)
                                      G.store[new.it]       <- G.non
      } 
    }
    close(pb)
    Gmax            <- seq_len(max(as.numeric(z.store)))
    returns         <- list(mu       = if(sw["mu.sw"])         mu.store[,Gmax,, drop=FALSE],
                            eta      = if(all(sw["s.sw"], Q0)) eta.store, 
                            load     = if(all(sw["l.sw"], Q0)) load.store[,,Gmax,, drop=FALSE], 
                            psi      = if(sw["psi.sw"])        psi.store[,Gmax,, drop=FALSE],
                            pi.prop  = if(sw["pi.sw"])         pi.store[Gmax,, drop=FALSE],
                            alpha    = if(not.fixed)           alpha.store,
                            discount = if(learn.d)             d.store,
                            rate     = if(MH.step)             mean(rate),
                            lab.rate = if(DP.lab.sw)           setNames(rowmeans(lab.rate), c("Move1", "Move2")),
                            z.store  = z.store,
                            ll.store = ll.store,
                            G.store  = G.store,
                            time     = init.time)
    return(returns)
  }