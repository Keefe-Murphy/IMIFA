################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MFA      <- function(Q, data, iters, N, P, G, 
                             sigma.mu, sigma.l, burnin, 
                             thinning, psi.alpha, psi.beta,
                             sw, verbose, clust, mu0g, ...) {
        
  # Define & initialise variables
    n.iters      <- round(max(iters), -1)
    n.store      <- length(iters)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    Gseq         <- seq_len(G)
    facnames     <- paste0("Factor ", seq_len(Q))
    gnames       <- paste0("Group ", Gseq)
    iternames    <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store   <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store  <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store   <- matrix(0, nr=G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store      <- matrix(0, nr=N, nc=n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    ll.store     <- rep(0, n.store)
    
    mu.sigma     <- 1/sigma.mu
    l.sigma      <- 1/sigma.l 
    pi.alpha     <- clust$pi.alpha
    z            <- clust$z
    pi.prop      <- t(prop.table(tabulate(z, nbins=G)))
    mu           <- do.call(cbind, lapply(Gseq, function(g) if(pi.prop[,g] > 0) colMeans(data[z == g,, drop=F]) else rep(0, P)))
    f            <- sim.f.p(N=N, Q=Q)
    lmat         <- lapply(Gseq, function(g) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F))
    psi.inv      <- do.call(cbind, lapply(Gseq, function(g) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    mu.zero      <- if(mu0g) mu else do.call(cbind, lapply(Gseq, function(g) colMeans(data)))
    if(round(sum(mu.zero)) == 0) {
      mu.zero    <- matrix(0, nr=1, nc=G)
    }
    if(Q > 0) {
      for(g in Gseq) {
        fact     <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } else                warning(paste0("Parameters of group ", g, " initialised by simulation from priors, not factanal: G=", G, ", Q=", Q), call.=F)
      }
    } else {
      psi.inv    <- do.call(cbind, lapply(Gseq, function(g) if(pi.prop[,g] > 0) 1/apply(data[z == g,, drop=F], 2, var) else rep(1, P)))
    }
    l.sigma      <- l.sigma * diag(Q)
    Qs           <- rep(Q, G)
    if(burnin     < 1)     {
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
      if(verbose) {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      nn         <- tabulate(z, nbins=G)
      z.ind      <- lapply(Gseq, function(g) z == g)
      
    # Means
      sum.data   <- lapply(Gseq, function(g) colSums(data[z.ind[[g]],, drop=F]))
      sum.f      <- lapply(Gseq, function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu         <- do.call(cbind, lapply(Gseq, function(g) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], 
                            P=P, sum.data=sum.data[[g]], sum.f=sum.f[[g]], lmat=lmat[[g]], mu.zero=mu.zero[,g])))
    
    # Scores & Loadings
      c.data     <- lapply(Gseq, function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(Q > 0)   {
        f        <- do.call(rbind, lapply(Gseq, function(g) sim.score(N=nn[g], lmat=lmat[[g]], 
                            c.data=c.data[[g]], psi.inv=psi.inv[,g], Q=Qs[g])))[obsnames,, drop=F]
        FtF      <- lapply(Gseq, function(g) crossprod(f[z.ind[[g]],, drop=F]))
        lmat     <- lapply(Gseq, function(g) sim.load(l.sigma=l.sigma, Q=Qs[g], c.data=c.data[[g]], 
                           P=P, f=f[z.ind[[g]],, drop=F], psi.inv=psi.inv[,g], FtF=FtF[[g]], shrink=F))
      }
                  
    # Uniquenesses
      psi.inv    <- do.call(cbind, lapply(Gseq, function(g) sim.psi.i(N=nn[g], P=P, psi.alpha=psi.alpha, 
                            psi.beta=psi.beta, c.data=c.data[[g]], f=f[z.ind[[g]],,drop=F], lmat=lmat[[g]])))
    
    # Mixing Proportions
      pi.prop    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi        <- 1/psi.inv
      Sigma      <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      z.res      <- sim.z(data=data, mu=mu, Sigma=Sigma, G=G, pi.prop=pi.prop)
      z          <- z.res$z
      
      if(is.element(iter, iters))  {
        new.it   <- which(iters == iter)
        log.like <- sum(z.res$log.likes) 
        if(sw["mu.sw"])             mu.store[,,new.it]      <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,new.it]       <- f
        if(all(sw["l.sw"], Q > 0)) {
          for(g in Gseq)  {
                                    load.store[,,g,new.it]  <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])            psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])             pi.store[,new.it]       <- pi.prop
                                    z.store[,new.it]        <- z 
                                    ll.store[new.it]        <- log.like
      }  
    }
    returns      <- list(mu       = if(sw["mu.sw"])            mu.store,
                         f        = if(all(sw["f.sw"], Q > 0)) f.store, 
                         load     = if(all(sw["l.sw"], Q > 0)) load.store, 
                         psi      = if(sw["psi.sw"])           psi.store,
                         pi.prop  = if(sw["pi.sw"])            pi.store,
                         z        = z.store,
                         ll.store = ll.store)
    attr(returns, "K")     <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
    return(returns)
  }