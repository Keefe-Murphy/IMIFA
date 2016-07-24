#######################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Dirichlet Process) ####
#######################################################################
  
# Gibbs Sampler Function
  gibbs.IMIFA       <- function(Q, data, iters, N, P, G, mu.zero, pp, sigma.l, MH.step,
                                sigma.mu, burnin, thinning, mu, trunc.G, MH.lower,
                                psi.alpha, psi.beta, verbose, gen.slice, alpha.d1,
                                alpha.dk, sw, cluster, phi.nu, b0, b1, prop, MH.upper,
                                beta.d1, beta.dk, adapt, epsilon, ...) {
        
  # Define & initialise variables
    n.iters         <- round(max(iters), -1)
    n.store         <- length(iters)
    Gs              <- seq_len(G)
    Ts              <- seq_len(trunc.G)
    Ps              <- seq_len(P)
    Ns              <- seq_len(N)
    obsnames        <- rownames(data)
    varnames        <- colnames(data)
    facnames        <- paste0("Factor ", seq_len(Q))
    gnames          <- paste0("Group ", Ts)
    iternames       <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store      <- array(0, dim=c(P, trunc.G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store       <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, facnames, iternames)
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
    non.empty       <- list()
    dimnames(z.store)      <- list(obsnames, iternames)
    dimnames(Q.store)      <- list(gnames, iternames)
    if(MH.step)    {
      rate          <- rep(0, n.store)
      alpha.store   <- rep(0, n.store)
    }
    
    mu.sigma        <- 1/sigma.mu
    z               <- cluster$z
    pi.alpha        <- cluster$pi.alpha
    pi.prop         <- cbind(cluster$pi.prop, sim.pi(pi.alpha=pi.alpha, nn=rep(0, trunc.G, inf.G=T))[,-Gs, drop=F])
    nn              <- tabulate(z, nbins=trunc.G)
    mu              <- cbind(mu, do.call(cbind, lapply(seq_len(trunc.G - G), function(g) sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero))))
    f               <- sim.f.p(N=N, Q=Q)
    phi             <- lapply(Ts, function(t) sim.phi.p(Q=Q, P=P, phi.nu=phi.nu))
    delta           <- lapply(Ts, function(t) c(sim.delta.p(alpha=alpha.d1, beta=beta.d1), sim.delta.p(Q=Q, alpha=alpha.dk, beta=beta.dk)))
    tau             <- lapply(delta, cumprod)
    lmat            <- lapply(Ts, function(t) matrix(unlist(lapply(Ps, function(j) sim.load.p(Q=Q, phi=phi[[t]][j,], tau=tau[[t]], P=P)), use.names=F), nr=P, byrow=T))
    psi.inv         <- do.call(cbind, lapply(Ts, function(t) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    for(g in which(nn > 2.5 * Q))      {
      fact          <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
      if(!inherits(fact, "try-error")) {
        f[z == g,]         <- fact$scores
        lmat[[g]]          <- fact$loadings
        psi.inv[,g]        <- 1/fact$uniquenesses
      } 
    }
    index           <- order(pi.prop, decreasing=TRUE)
    pi.prop         <- pi.prop[,index, drop=F]
    mu              <- mu[,index, drop=FALSE]
    lmat            <- lmat[index]
    psi.inv         <- psi.inv[,index, drop=FALSE]
    csi             <- pp * (1 - pp)^(Ts - 1)
    if(burnin        < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- array(unlist(lmat, use.names=F), dim=c(P, Q, G))
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gs, N=N, pi.prop=pi.prop, Sigma=lapply(Gs,
                                  function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])), Q0=Qs > 0)$log.likes)
      Q.store[,1]          <- Qs
      G.store[1]           <- G
      if(MH.step)  {
        rate[1]            <- 0
        alpha.store[1]     <- pi.alpha 
      }
    }
    
  # Iterate
    for(iter in seq_len(max(iters))[-1]) { 
      if(verbose)  {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
    
    # Slice Sampler
      csi           <- if(gen.slice) csi else pi.prop
      u.slice       <- runif(N, 0, csi[z])
      Gs            <- seq_len(max(unlist(lapply(Ns, function(i) sum(u.slice[i] < pi.prop)))))
      slice.ind     <- do.call(cbind, lapply(Gs, function(g, x=csi[g]) (u.slice < x)/x))
    
    # Mixing Proportions
      weights       <- sim.pi(pi.alpha=pi.alpha, nn=nn, inf.G=T, len=trunc.G)
      pi.prop       <- weights$pi.prop
      if(MH.step)  {
        Vs          <- weights$Vs
      }
    
    # Cluster Labels
      psi           <- 1/psi.inv
      Sigma         <- lapply(Gs, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      Q0            <- Qs > 0
      Q1            <- Qs > 1
      z.res         <- sim.z(data=data, mu=mu[,Gs], Sigma=Sigma, Gseq=Gs, N=N, pi.prop=pi.prop[,Gs, drop=F], slice.ind=slice.ind, Q0=Q0[Gs])
      z             <- z.res$z
      nn            <- tabulate(z, nbins=trunc.G)
      nn0           <- nn > 0
      z.ind         <- lapply(Gs, function(g) z == g)
      dat.g         <- lapply(Gs, function(g) data[z.ind[[g]],, drop=F])
      
    # Scores & Loadings
      c.data        <- lapply(Gs, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        f          <- matrix(, nr=N, nc=0)
        f.tmp      <- lapply(Gs, function(g) f[z.ind[[g]],, drop=F])
        lmat       <- lapply(Gs, function(g) matrix(, nr=P, nc=0))
      } else {
        f.tmp      <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) sim.score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else matrix(, nr=ifelse(Q0[g], 0, nn[g]), nc=Qs[g]))
        FtF        <- lapply(Gs, function(g) if(nn0[g]) crossprod(f.tmp[[g]]))
        lmat[Gs]   <- lapply(Gs, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Ps, function(j) sim.load(Q=Qs[g], P=P, c.data=c.data[[g]][,j], FtF=FtF[[g]],
                             f=f.tmp[[g]], psi.inv=psi.inv[,g][j], Q1=Q1[g], phi=phi[[g]][j,], tau=tau[[g]], shrink=T)), use.names=F), nr=P, byrow=T) else 
                             matrix(unlist(lapply(Ps, function(j) sim.load.p(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], P=P)), use.names=F), nr=P, byrow=F))
        f.tmp      <- if(length(unique(Qs)) != 1) lapply(Gs, function(g) cbind(f.tmp[[g]], matrix(0, nr=nn[g], nc=max(Qs) - Qs[g]))) else f.tmp
        q0ng       <- !Q0 & nn0
        if(any(q0ng)) {
          f.tmp[q0ng]      <- lapply(Gs[q0ng], function(g, x=f.tmp[[g]]) { row.names(x) <- obsnames[z.ind[[g]]]; x })
        }
        f          <- do.call(rbind, f.tmp)[obsnames,, drop=F]
      }
      
    # Means
      sum.data      <- lapply(dat.g, colSums)
      sum.f         <- lapply(f.tmp, colSums)
      mu[,Gs]       <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.data=sum.data[[g]], 
                               sum.f=sum.f[[g]][seq_len(Qs[g])], lmat=lmat[[g]], mu.zero=mu.zero) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero)))
      
    # Uniquenesses
      psi.inv[,Gs]  <- do.call(cbind, lapply(Gs, function(g) if(nn0[g]) sim.psi.i(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta, 
                               P=P, f=f.tmp[[g]][,seq_len(Qs[g]), drop=F], lmat=lmat[[g]]) else sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    
    # Local Shrinkage
      load.2        <- lapply(lmat[Gs], function(lg) lg * lg)
      phi[Gs]       <- lapply(Gs, function(g) if(nn0[g]) sim.phi(Q=Qs[g], P=P, phi.nu=phi.nu, 
                       tau=tau[[g]], load.2=load.2[[g]]) else sim.phi.p(Q=Qs[g], P=P, phi.nu=phi.nu))
    
    # Global Shrinkage
      sum.terms    <- lapply(Gs, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gs)   {
        Qg         <- Qs[g]
        nn0g       <- nn0[g]
        if(nn0g)   {
          sumtermg <- sum.terms[[g]]  
        }
      if(Q0[g])  {
        delta[[g]][1]    <- if(nn0g) sim.delta1(Q=Qg, alpha.d1=alpha.d1, delta.1=delta[[g]][1], P=P, beta.d1=beta.d1, 
                            tau=tau[[g]], sum.term=sumtermg) else sim.delta.p(alpha=alpha.d1, beta=beta.d1)
        tau[[g]]         <- cumprod(delta[[g]])
      }
      if(Q1[g])  {
        for(k in seq_len(Qg)[-1]) { 
          delta[[g]][k]  <- if(nn0g) sim.deltak(Q=Qg, alpha.dk=alpha.dk, delta.k=delta[[g]][k], P=P, beta.dk=beta.dk, k=k, 
                            tau.kq=tau[[g]][k:Qg], sum.term.kq=sumtermg[k:Qg]) else sim.delta.p(alpha=alpha.dk, beta=beta.dk)
          tau[[g]]       <- cumprod(delta[[g]])
        }
      }
    }
    
    # Adaptation  
      if(all(adapt, iter > burnin)) {      
        if(runif(1) < 1/exp(b0 + b1 * (iter - burnin))) {
          lind      <- lapply(Gs, function(g) if(Q0[g]) colSums(abs(lmat[[g]]) < epsilon)/P else 0)
          colvec    <- lapply(lind, function(lx) lx >= prop)
          nonred    <- lapply(colvec, function(cv) which(cv == 0))
          numred    <- lapply(colvec, sum)
          notred    <- unlist(lapply(Gs, function(g) numred[[g]] == 0 && nn0[g]), use.names=F)
          Qs.old    <- Qs
          Qs[Gs]    <- unlist(lapply(Gs, function(g) if(notred[g]) Qs.old[g] + 1 else Qs.old[g] - numred[[g]]), use.names=F)
          phi[Gs]   <- lapply(Gs, function(g) if(notred[g]) cbind(phi[[g]][,seq_len(Qs.old[g])], rgamma(n=P, shape=phi.nu, rate=phi.nu)) else phi[[g]][,nonred[[g]], drop=F])
          delta[Gs] <- lapply(Gs, function(g) if(notred[g]) c(delta[[g]][seq_len(Qs.old[g])], rgamma(n=1, shape=alpha.dk, rate=beta.dk)) else delta[[g]][nonred[[g]]])  
          tau[Gs]   <- lapply(delta[Gs], cumprod)
          lmat[Gs]  <- lapply(Gs, function(g) if(notred[g]) cbind(lmat[[g]][,seq_len(Qs.old[g])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[g]], drop=F])
          f         <- if(max(Qs) > max(Qs.old)) cbind(f[,seq_len(max(Qs.old))], rnorm(N)) else f[,seq_len(max(Qs)), drop=F]
          Q.max     <- max(Qs[nn0])
          if(max(Qs[!nn0], 0) > Q.max) {
            Qs[Qs > Q.max] <- Q.max
            Q.maxseq       <- seq_len(max(Qs))
            for(t0 in Ts[!nn0])        {  
              phi[[t0]]    <- phi[[t0]][,Q.maxseq, drop=F]
              delta[[t0]]  <- delta[[t0]][Q.maxseq, drop=F]
              tau[[t0]]    <- tau[[t0]][Q.maxseq, drop=F]
              lmat[[t0]]   <- lmat[[t0]][,Q.maxseq, drop=F]
            }
            f       <- f[,Q.maxseq, drop=F]
          }
        }
      }
    
    # Alpha
      if(MH.step)   {
        MH.alpha    <- sim.alpha(lower=MH.lower, upper=MH.upper, trunc.G=trunc.G, alpha=pi.alpha, Vs=Vs) 
        pi.alpha    <- MH.alpha$alpha
      }
    
    if(any(Qs > Q.star))      stop(paste0("Q cannot exceed initial number of loadings columns: try increasing range.Q from ", Q.star))
      if(is.element(iter, iters))   {
        new.it      <- which(iters == iter)
        log.like    <- sum(z.res$log.likes)
        if(sw["mu.sw"])    mu.store[,,new.it]       <- mu 
        if(all(sw["f.sw"], 
           any(Q0)))  f.store[,seq_len(max(Qs)),new.it]    <- f
        if(sw["l.sw"])   {
          for(g in Gs)   {
            if(Q0[g]) load.store[,seq_len(Qs[g]),g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"]) psi.store[,,new.it]        <- psi
        if(sw["pi.sw"])  pi.store[,new.it]          <- pi.prop
        if(MH.step) {    rate[new.it]               <- MH.alpha$rate
                         alpha.store[new.it]        <- pi.alpha }
                         z.store[,new.it]           <- z 
                         ll.store[new.it]           <- log.like
                         Q.store[,new.it]           <- Qs
                         G.store[new.it]            <- sum(nn0)
                         non.empty[[new.it]]        <- which(nn0)
      } 
    }
  
    returns         <- list(mu       = if(sw["mu.sw"])  mu.store,
                            f        = if(sw["f.sw"])   as.simple_sparse_array(f.store), 
                            load     = if(sw["l.sw"])   as.simple_sparse_array(load.store), 
                            psi      = if(sw["psi.sw"]) psi.store,
                            pi.prop  = if(sw["pi.sw"])  pi.store,
                            rate     = if(MH.step)      mean(rate),
                            alpha    = if(MH.step)      alpha.store,
                            z.store  = z.store,
                            ll.store = ll.store,
                            Q.store  = Q.store,
                            G.store  = G.store,
                            nonempty = non.empty)
    return(returns)
  }