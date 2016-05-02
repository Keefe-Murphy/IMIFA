################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MIFA      <- function(Q, data, iters, N, P, G, 
                              sigma.mu, burnin, thinning, 
                              psi.alpha, psi.beta, verbose, 
                              sw, clust, mu0g, phi.nu, b0,
                              b1, alpha.d1, alpha.dk, beta.d1,
                              beta.dk, adapt, prop, epsilon, ...) {
        
  # Define & initialise variables
    n.iters       <- round(max(iters), -1)
    n.store       <- length(iters)
    obsnames      <- rownames(data)
    varnames      <- colnames(data)
    Gseq          <- seq_len(G)
    facnames      <- paste0("Factor ", seq_len(Q))
    gnames        <- paste0("Group ", Gseq)
    iternames     <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store    <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store     <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store  <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store   <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store    <- matrix(0, nr=G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store       <- matrix(0, nr=N, nc=n.store)
    ll.store      <- rep(0, n.store)
    Q.star        <- Q
    Qs            <- rep(Q, G)
    Q.store       <- matrix(0, nr=G, nc=n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    dimnames(Q.store)      <- list(gnames, iternames)
    
    mu.sigma      <- 1/sigma.mu
    l.sigma       <- 1/sigma.l 
    pi.alpha      <- clust$pi.alpha
    z             <- clust$z
    pi.prop       <- t(prop.table(tabulate(z, nbins=G)))
    mu            <- do.call(cbind, lapply(Gseq, function(g) if(pi.prop[,g] > 0) colMeans(data[z == g,, drop=F]) else rep(0, P)))
    f             <- sim.f.p(N=N, Q=Q)
    phi           <- lapply(Gseq, function(g) sim.phi.p(Q=Q, P=P, phi.nu=phi.nu))
    delta         <- lapply(Gseq, function(g) sim.delta.p(Q=Q, alpha.d1=alpha.d1, alpha.dk=alpha.dk, beta.d1=beta.d1, beta.dk=beta.dk))
    tau           <- lapply(delta, cumprod)
    lmat          <- lapply(Gseq, function(g) sim.load.p(Q=Q, phi=phi[[g]], tau=tau[[g]], P=P))
    psi.inv       <- do.call(cbind, lapply(Gseq, function(g) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    mu.zero       <- if(mu0g) mu else do.call(cbind, lapply(Gseq, function(g) colMeans(data)))
    if(round(sum(mu.zero)) == 0) {
      mu.zero     <- matrix(0, nr=1, nc=G)
    }
    for(g in Gseq) {
      fact        <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
      if(!inherits(fact, "try-error")) {
        f[z == g,]         <- fact$scores
        lmat[[g]]          <- fact$loadings
        psi.inv[,g]        <- 1/fact$uniquenesses
      } else                  warning(paste0("Parameters of group ", g, " initialised by simulation from priors, not factanal: G=", G, ", Q=", Q), call.=F)
    }
    if(burnin      < 1)     {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      for(g in Gseq) {
        load.store[,,g,1]  <- lmat[[g]]
      }
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, G=G, pi.prop=pi.prop, Sigma=lapply(Gseq,
                                  function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])))$log.likes)
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
      nn          <- tabulate(z, nbins=G)
      z.ind       <- lapply(Gseq, function(g) z == g)
      
    # Means
      sum.data    <- lapply(Gseq, function(g) colSums(data[z.ind[[g]],, drop=F]))
      sum.f       <- lapply(Gseq, function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu          <- do.call(cbind, lapply(Gseq, function(g) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, 
                             sum.data=sum.data[[g]], sum.f=sum.f[[g]][seq_len(Qs[g])], lmat=lmat[[g]], mu.zero=mu.zero[,g])))
    
    # Scores & Loadings
      c.data      <- lapply(Gseq, function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(all(Qs   == 0))  {
        f         <- matrix(, nr=N, nc=0)
        lmat      <- lapply(Gseq, function(g) matrix(, nr=P, nc=0))
      } else {
        fg        <- lapply(Gseq, function(g) matrix(0, nr=nn[g], nc=max(Qs)))
        fgnames   <- lapply(Gseq, function(g) obsnames[z.ind[[g]]])
        fg        <- mapply(function(mats, nams) {rownames(mats) <- nams; mats}, fg, fgnames)
        for(g in Gseq)    {
          Qg      <- Qs[g]
          Qgs     <- seq_len(Qg)
          nng     <- nn[g]
          c.datg  <- c.data[[g]]
          psi.ig  <- psi.inv[,g]
          if(Qg    > 0)   {
            fgg            <- sim.score(N=nng, lmat=lmat[[g]], Q=Qg, 
                                        c.data=c.datg, psi.inv=psi.ig)
            FtF            <- crossprod(fgg)
            lmat[[g]]      <- sim.load(Q=Qg, c.data=c.datg, P=P, f=fgg, FtF=FtF,
                                       psi.inv=psi.ig, phi=phi[[g]], tau=tau[[g]])
            fg[[g]][,Qgs]  <- fgg
          } else {
            fg[[g]][,Qgs]  <- matrix(, nr=nng, nc=0)
            lmat[[g]]      <- matrix(, nr=P, nc=0)
          }
        }
        f         <- do.call(rbind, fg)[obsnames,, drop=F]
      }
                  
    # Uniquenesses
      psi.inv     <- do.call(cbind, lapply(Gseq, function(g) sim.psi.i(N=nn[g], psi.alpha=psi.alpha, psi.beta=psi.beta, 
                             P=P, c.data=c.data[[g]], f=f[z.ind[[g]],seq_len(Qs[g]),drop=F], lmat=lmat[[g]])))
    
    # Local Shrinkage
      load.2      <- lapply(lmat, function(lg) lg * lg)
      phi         <- lapply(Gseq, function(g) sim.phi(Q=Qs[g], P=P, phi.nu=phi.nu, tau=tau[[g]], load.2=load.2[[g]]))
    
    # Global Shrinkage
      sum.terms   <- lapply(Gseq, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gseq)    {
        Qg        <- Qs[g]
        sum.termg <- sum.terms[[g]]
        if(Qg      > 0) {
          delta[[g]][1]    <- sim.delta1(Q=Qg, alpha.d1=alpha.d1, delta=delta[[g]], P=P,
                                         beta.d1=beta.d1, tau=tau[[g]], sum.term=sum.termg)
          tau[[g]]         <- cumprod(delta[[g]])
        }
        if(Qg      > 1) {
          for(k in seq_len(Qg)[-1]) { 
            delta[[g]][k]  <- sim.deltak(Q=Qg, alpha.dk=alpha.dk, delta=delta[[g]], P=P,
                                         beta.dk=beta.dk, k=k, tau=tau[[g]], sum.term=sum.termg)
            tau[[g]]       <- cumprod(delta[[g]])
          }
        }
      }
    
    # Adaptation  
      if(all(adapt, iter > burnin)) {      
        prob      <- 1/exp(b0 + b1 * pmax(iter - burnin, 0))
        unif      <- runif(n=1, min=0, max=1)     
        if(unif    < prob) { 
          lind    <- lapply(Gseq, function(g) if(Qs[g] > 0) colSums(abs(lmat[[g]]) < epsilon)/P else 0)
          colvec  <- lapply(lind, function(lx) lx >= prop)
          nonred  <- lapply(colvec, function(cv) which(cv == 0))
          numred  <- lapply(colvec, sum)
          notred  <- unlist(lapply(Gseq, function(g) numred[[g]] == 0))
          Qs.old  <- Qs
          Qs      <- unlist(lapply(Gseq, function(g) if(notred[g]) Qs.old[g] + 1 else Qs.old[g] - numred[[g]]))
          phi     <- lapply(Gseq, function(g) if(notred[g]) cbind(phi[[g]][,seq_len(Qs.old[g])], rgamma(n=P, shape=phi.nu, rate=phi.nu)) else phi[[g]][,nonred[[g]], drop=F])
          delta   <- lapply(Gseq, function(g) if(notred[g]) c(delta[[g]][seq_len(Qs.old[g])], rgamma(n=1, shape=alpha.dk, rate=beta.dk)) else delta[[g]][nonred[[g]]])  
          tau     <- lapply(delta, cumprod)
          lmat    <- lapply(Gseq, function(g) if(notred[g]) cbind(lmat[[g]][,seq_len(Qs.old[g])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qs[g]] * tau[[g]][Qs[g]])))) else lmat[[g]][,nonred[[g]], drop=F])
          f       <- if(max(Qs) > max(Qs.old)) cbind(f[,seq_len(max(Qs.old))], rnorm(n=N, mean=0, sd=1)) else f[,seq_len(max(Qs)), drop=F]
        }
      }
    
    # Mixing Proportions
      pi.prop     <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi         <- 1/psi.inv
      Sigma       <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      z.res       <- sim.z(data=data, mu=mu, Sigma=Sigma, G=G, pi.prop=pi.prop)
      z           <- z.res$z
      
    if(any(Qs > Q.star))      stop(paste0("Q cannot exceed initial number of loadings columns: try increasing Q.star from ", Q.star))
      if(is.element(iter, iters))  {
        new.it    <- which(iters == iter)
        log.like  <- sum(z.res$log.likes)
        if(sw["mu.sw"])       mu.store[,,new.it]   <- mu  
        if(all(sw["f.sw"], 
           any(Qs > 0)))      f.store[,seq_len(max(Qs)),new.it]    <- f
        if(sw["l.sw"])   {
          for(g in Gseq) {
            if(Qs[g] > 0)     load.store[,seq_len(Qs[g]),g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])      psi.store[,,new.it]  <- psi
        if(sw["pi.sw"])       pi.store[,new.it]    <- pi.prop
                              z.store[,new.it]     <- z 
                              ll.store[new.it]     <- log.like  
                              Q.store[,new.it]     <- Qs
      }
    }
    returns       <- list(mu       = if(sw["mu.sw"])  mu.store,
                          f        = if(sw["f.sw"])   as.simple_sparse_array(f.store), 
                          load     = if(sw["l.sw"])   as.simple_sparse_array(load.store), 
                          psi      = if(sw["psi.sw"]) psi.store,
                          pi.prop  = if(sw["pi.sw"])  pi.store,
                          z        = z.store,
                          ll.store = ll.store,
                          Q.store  = Q.store)
    return(returns)
  }