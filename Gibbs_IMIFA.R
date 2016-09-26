#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  gibbs.IMIFA       <- function(Q, data, iters, N, P, G, mu.zero, rho, sigma.l, alpha.step, mu, sw, 
                                sigma.mu, burnin, thinning, trunc.G, a.hyper, psi.alpha, psi.beta, 
                                verbose, gen.slice, alpha.d1, discount, alpha.dk, cluster, b0, b1, adapt,
                                phi.nu, prop, d.hyper, beta.d1, beta.dk, adapt.at, epsilon, learn.d, ...) {
        
  # Define & initialise variables
    start.time      <- proc.time()
    total           <- max(iters)
    if(verbose)        pb  <- txtProgressBar(min=0, max=total, style=3)
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
    if(sw["mu.sw"])  {
      mu.store      <- array(0, dim=c(P, trunc.G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["s.sw"])   {
      eta.store     <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)  <- list(obsnames, facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store    <- array(0, dim=c(P, Q, trunc.G, n.store))
      dimnames(load.store) <- list(varnames, facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store     <- array(0, dim=c(P, trunc.G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store      <- matrix(0, nr=trunc.G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store         <- matrix(0, nr=N, nc=n.store)
    ll.store        <- rep(0, n.store)
    Q.star          <- Q
    Qs              <- rep(Q, trunc.G)
    Q.store         <- matrix(0, nr=trunc.G, nc=n.store)
    G.store         <- rep(0, n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    dimnames(Q.store)      <- list(gnames, iternames)
    not.fixed       <- alpha.step != "fixed"
    if(not.fixed) {
      alpha.store   <- rep(0, n.store)
      alpha.shape   <- a.hyper[1]
      alpha.rate    <- a.hyper[2]
    }
    if(learn.d)   {
      d.store       <- rep(0, n.store)
      d.shape1      <- d.hyper[1]
      d.shape2      <- d.hyper[2]
    }
    MH.step         <- alpha.step == "metropolis"
    if(MH.step)   {
      rate          <- rep(0, n.store)
    }
    mu.sigma        <- 1/sigma.mu
    z               <- cluster$z
    pi.alpha        <- cluster$pi.alpha
    pi.prop         <- c(cluster$pi.prop, sim.pi(pi.alpha=pi.alpha, nn=rep(0, trunc.G), N=N, inf.G=TRUE, len=trunc.G, discount=discount)$pi.prop[-Gs])
    nn              <- tabulate(z, nbins=trunc.G)
    mu              <- cbind(mu, vapply(seq_len(trunc.G - G), function(g) sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P)))
    eta             <- sim.eta.p(N=N, Q=Q)
    phi             <- lapply(Ts, function(t) sim.phi.p(Q=Q, P=P, phi.nu=phi.nu))
    delta           <- lapply(Ts, function(t) c(sim.delta.p(alpha=alpha.d1, beta=beta.d1), sim.delta.p(Q=Q, alpha=alpha.dk, beta=beta.dk)))
    tau             <- lapply(delta, cumprod)
    lmat            <- lapply(Ts, function(t) matrix(unlist(lapply(Ps, function(j) sim.load.p(Q=Q, phi=phi[[t]][j,], tau=tau[[t]], P=P)), use.names=FALSE), nr=P, byrow=TRUE))
    psi.inv         <- vapply(Ts, function(t) sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    for(g in which(nn > 2.5 * Q))      {
      fact          <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
      if(!inherits(fact, "try-error")) {
        eta[z == g,]       <- fact$scores
        lmat[[g]]          <- fact$loadings
        psi.inv[,g]        <- 1/fact$uniquenesses
      } 
    }
    index           <- order(pi.prop, decreasing=TRUE)
    pi.prop         <- pi.prop[index]
    mu              <- mu[,index, drop=FALSE]
    phi             <- phi[index]
    delta           <- delta[index]
    tau             <- tau[index]
    lmat            <- lmat[index]
    psi.inv         <- psi.inv[,index, drop=FALSE]
    ksi             <- (1 - rho) * rho^(Ts - 1)
    k.x             <- .Machine$double.xmin
    ksi[ksi < k.x]  <- k.x
    if(burnin        < 1)  {
      mu.store[,,1]        <- mu
      eta.store[,,1]       <- eta
      load.store[,,,1]     <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, trunc.G))
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop/sum(pi.prop)
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu[,Gs], Gseq=Gs, N=N, pi.prop=pi.store[Gs,1], sigma=lapply(Gs,
                                  function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])), Q0=Qs > 0)$log.likes)
      Q.store[,1]          <- Qs
      G.store[1]           <- sum(nn > 0)
      if(not.fixed) {
        alpha.store[1]     <- pi.alpha 
      }
      if(learn.d)   {
        d.store[1]         <- discount
      }
    }
    init.time       <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose    && iter < burnin) setTxtProgressBar(pb, iter)
    
    # Mixing Proportions
      weights       <- sim.pi(pi.alpha=pi.alpha, nn=nn, N=N, inf.G=TRUE, len=trunc.G, discount=discount)
      pi.prop       <- weights$pi.prop
      if(MH.step)      Vs  <- weights$Vs
      
    # Slice Sampler
      if(!gen.slice) {
        index       <- order(pi.prop, decreasing=TRUE)
        pi.prop     <- ksi <- pi.prop[index]
        mu          <- mu[,index, drop=FALSE]
        phi         <- phi[index]
        delta       <- delta[index]
        tau         <- tau[index]
        lmat        <- lmat[index]
        Qs          <- Qs[index]
        psi.inv     <- psi.inv[,index, drop=FALSE]
        ksi[ksi < k.x]     <- k.x
      }
      u.slice       <- runif(N, 0, ksi[z])
      Gs            <- seq_len(max(vapply(Ns, function(i) sum(u.slice[i] < pi.prop), numeric(1))))
      slice.ind     <- vapply(Gs, function(g, x=ksi[g]) (u.slice < x)/x, numeric(N))
    
    # Cluster Labels
      psi           <- 1/psi.inv
      sigma         <- lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      Q0            <- Qs[Gs] > 0
      Q1            <- Qs[Gs] > 1
      z.res         <- sim.z(data=data, mu=mu[,Gs, drop=FALSE], sigma=sigma, Gseq=Gs, N=N, pi.prop=pi.prop[Gs], slice.ind=slice.ind, Q0=Q0[Gs])
      z             <- z.res$z
      nn            <- tabulate(z, nbins=trunc.G)
      nn0           <- nn > 0
      z.ind         <- lapply(Gs, function(g) z == g)
      dat.g         <- lapply(Gs, function(g) data[z.ind[[g]],, drop=FALSE])
      
    # Scores & Loadings
      c.data        <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta         <- matrix(, nr=N, nc=0)
        eta.tmp     <- lapply(Gs, function(g) eta[z.ind[[g]],, drop=FALSE])
        lmat[Gs]    <- lapply(Gs, function(g) matrix(, nr=P, nc=0))
      } else {
        eta.tmp     <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) sim.score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else matrix(, nr=ifelse(Q0[g], 0, nn[g]), nc=Qs[g]))
        EtE         <- lapply(Gs, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat[Gs]    <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Ps, function(j) sim.load(Q=Qs[g], P=P, c.data=c.data[[g]][,j], EtE=EtE[[g]],
                             eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], Q1=Q1[g], phi=phi[[g]][j,], tau=tau[[g]], shrink=TRUE)), use.names=FALSE), nr=P, byrow=TRUE) else 
                             matrix(unlist(lapply(Ps, function(j) sim.load.p(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], P=P)), use.names=FALSE), nr=P, byrow=FALSE))
        eta.tmp     <- if(length(unique(Qs)) != 1) lapply(Gs, function(g) cbind(eta.tmp[[g]], matrix(0, nr=nn[g], nc=max(Qs) - Qs[g]))) else eta.tmp
        q0ng        <- !Q0 & nn0[Gs]
        if(any(q0ng)) {
          eta.tmp[q0ng]    <- lapply(Gs[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- obsnames[z.ind[[g]]]; x })
        }
        eta         <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      }
      
    # Means
      sum.data      <- lapply(dat.g, colSums)
      sum.eta       <- lapply(eta.tmp, colSums)
      mu[,Gs]       <- vapply(Gs, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])],
                              sum.data=sum.data[[g]], lmat=lmat[[g]], mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero), numeric(P))
      
    # Uniquenesses
      psi.inv[,Gs]  <- vapply(Gs, function(g) if(nn0[g]) sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, lmat=lmat[[g]],
                              P=P, eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta), numeric(P))
    
    # Local Shrinkage
      load.2        <- lapply(lmat[Gs], function(lg) lg * lg)
      phi[Gs]       <- lapply(Gs, function(g) if(nn0[g]) sim.phi(Q=Qs[g], P=P, phi.nu=phi.nu, 
                       tau=tau[[g]], load.2=load.2[[g]]) else sim.phi.p(Q=Qs[g], P=P, phi.nu=phi.nu))
    
    # Global Shrinkage
      sum.terms    <- lapply(Gs, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gs) {
        Qg         <- Qs[g]
        if(nn0[g]) {
          for(k in seq_len(Qg)) { 
            delta[[g]][k]  <- if(k > 1) sim.deltak(alpha.dk=alpha.dk, beta.dk=beta.dk, delta.k=delta[[g]][k], Q=Qg, P=P, 
                              k=k, tau.kq=tau[[g]][k:Qg], sum.term.kq=sum.terms[[g]][k:Qg]) else sim.delta1(alpha.d1=alpha.d1,
                              beta.d1=beta.d1, delta.1=delta[[g]][1], Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]])
            tau[[g]]       <- cumprod(delta[[g]])
          }
        } else {
          for(k in seq_len(Qg)) { 
            delta[[g]][k]  <- if(k > 1) sim.delta.p(alpha=alpha.dk, beta=beta.dk) else sim.delta.p(alpha=alpha.d1, beta=beta.d1)
            tau[[g]]       <- cumprod(delta[[g]])
          }
        }
      }
    
    # Adaptation  
      if(all(adapt, iter > adapt.at)) {      
        if(runif(1)  < ifelse(iter < burnin, 0.5, 1/exp(b0 + b1 * (iter - adapt.at)))) {
          lind      <- lapply(Gs, function(g) if(all(Q0[g], nn0[g])) colSums(abs(lmat[[g]]) < epsilon)/P else rep(0, Qs[g]))
          colvec    <- lapply(lind, function(lx) lx >= prop)
          nonred    <- lapply(colvec, function(cv) which(cv == 0))
          numred    <- lapply(colvec, sum)
          notred    <- vapply(Gs, function(g) numred[[g]] == 0, logical(1))
          Qs.old    <- Qs[Gs]
          Qs[Gs]    <- vapply(Gs, function(g) if(notred[g]) Qs.old[g] + 1 else Qs.old[g] - numred[[g]], numeric(1))
          phi[Gs]   <- lapply(Gs, function(g) if(notred[g]) cbind(phi[[g]][,seq_len(Qs.old[g])], rgamma(n=P, shape=phi.nu, rate=phi.nu)) else phi[[g]][,nonred[[g]], drop=FALSE])
          delta[Gs] <- lapply(Gs, function(g) if(notred[g]) c(delta[[g]][seq_len(Qs.old[g])], rgamma(n=1, shape=alpha.dk, rate=beta.dk)) else delta[[g]][nonred[[g]]])  
          tau[Gs]   <- lapply(delta[Gs], cumprod)
          lmat[Gs]  <- lapply(Gs, function(g, Qg=Qs[g]) if(notred[g]) cbind(lmat[[g]][,seq_len(Qs.old[g])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qg] * tau[[g]][Qg])))) else lmat[[g]][,nonred[[g]], drop=FALSE])
          Qemp      <- Qs[!nn0]
          Qmax      <- max(Qs[nn0])
          Qmaxseq   <- seq_len(Qmax)
          Qmaxold   <- max(Qs.old[nn0])
          eta       <- if(Qmax > Qmaxold) cbind(eta[,seq_len(Qmaxold)], rnorm(N)) else eta[,Qmaxseq, drop=FALSE]
          if(Qmax    < max(Qemp, 0)) {
            Qs[Qmax  < Qs] <- Qmax
            for(t  in  Ts[!nn0][Qemp > Qmax]) {  
              phi[[t]]     <- phi[[t]][,Qmaxseq,  drop=FALSE]
              delta[[t]]   <- delta[[t]][Qmaxseq]
              tau[[t]]     <- tau[[t]][Qmaxseq]
              lmat[[t]]    <- lmat[[t]][,Qmaxseq, drop=FALSE]
            }
          }
        }
      }
    
    # Alpha
      if(not.fixed) {
        if(MH.step) {
          MH.alpha  <- sim.alpha.m(alpha=pi.alpha, lower=alpha.shape, upper=alpha.rate, trunc.G=trunc.G, Vs=Vs, discount=discount) 
          pi.alpha  <- MH.alpha$alpha  
        } else {
          pi.alpha  <- sim.alpha.g(alpha=pi.alpha, shape=alpha.shape, rate=alpha.rate, G=G, N=N, discount=discount) 
        }
      }
    
    if(any(Qs > Q.star))      stop(paste0("Q cannot exceed initial number of loadings columns: try increasing range.Q from ", Q.star))
      if(is.element(iter, iters))   {
        if(verbose)    setTxtProgressBar(pb, iter)
        new.it      <- which(iters == iter)
        if(sw["mu.sw"])  mu.store[,,new.it]         <- mu  
        if(all(sw["s.sw"], 
           any(Q0)))   eta.store[,seq_len(max(Qs)),new.it]  <- eta 
        if(sw["l.sw"])   {
          for(g in Gs)   {
            if(Q0[g])  load.store[,seq_len(Qs[g]),g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"]) psi.store[,,new.it]        <- psi
        if(sw["pi.sw"])  pi.store[,new.it]          <- pi.prop
        if(not.fixed)    alpha.store[new.it]        <- pi.alpha
        if(learn.d)      d.store[new.it]            <- discount
        if(MH.step)      rate[new.it]               <- MH.alpha$rate
                         z.store[,new.it]           <- z 
                         ll.store[new.it]           <- sum(z.res$log.likes)
                         Q.store[,new.it]           <- Qs
                         G.store[new.it]            <- sum(nn0)
      } 
    }
    close(pb)
    Gmax            <- seq_len(max(as.numeric(z.store)))
    Qmax            <- seq_len(max(Q.store))
    returns         <- list(mu       = if(sw["mu.sw"])  mu.store[,Gmax,, drop=FALSE],
                            eta      = if(sw["s.sw"])   as.simple_sparse_array(eta.store[,Qmax,, drop=FALSE]), 
                            load     = if(sw["l.sw"])   as.simple_sparse_array(load.store[,Qmax,Gmax,, drop=FALSE]), 
                            psi      = if(sw["psi.sw"]) psi.store[,Gmax,, drop=FALSE],
                            pi.prop  = if(sw["pi.sw"])  pi.store[Gmax,, drop=FALSE],
                            alpha    = if(not.fixed)    alpha.store,
                            discount = if(learn.d)      discount,
                            rate     = if(MH.step)      mean(rate),
                            z.store  = z.store,
                            ll.store = ll.store,
                            G.store  = G.store,
                            Q.store  = Q.store[Gmax,, drop=FALSE],
                            time     = init.time)
    return(returns)
  }