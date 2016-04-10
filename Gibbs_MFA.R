################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MFA      <- function(Q = NULL, data = NULL, n.iters = NULL,
                             N = NULL, P = NULL, sigma.mu = NULL,
                             psi.alpha = NULL, psi.beta = NULL, G = NULL,
                             burnin = NULL, thinning = NULL, sw = NULL,
                             n.store = NULL, verbose = NULL, sigma.l = NULL,
                             alpha.pi = NULL, zinit = NULL, zlist = NULL,  ...) {
        
  # Define & initialise variables
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    gnames       <- paste0("Group ", seq_len(G))
    iternames    <- paste0("Iteration", seq_len(n.store))
    iters        <- seq(from=burnin + 1, to=n.iters, by=thinning)
    iters        <- iters[iters > 0]
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
    
    sigma.mu     <- 1/sigma.mu
    sigma.l      <- 1/sigma.l
    pi.alpha     <- rep(alpha.pi, G)
    if(zinit == "priors") {
      pi.prop    <- sim.pi(pi.alpha=pi.alpha)
      z          <- sim.z.p(N=N, prob.z=pi.prop)
    } else   {
      if(zinit == "list") {
        z        <- zlist
      } else {
        z        <- factor(kmeans(data, G, nstart=100)$cluster, levels=seq_len(G))
      }
      pi.prop    <- prop.table(tabulate(z, nbins=G))
    } 
    zinit        <- z
    mu           <- sim.mu.mp(P=P, sigma.mu=sigma.mu, G=G) 
    f            <- sim.f.mp(Q=Q, N=N)
    lmat         <- sim.load.mp(Q=Q, P=P, sigma.l=sigma.l, G=G)
    psi.inv      <- sim.psi.imp(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, G=G)
    l.sigma      <- sigma.l * diag(Q)
    Qs           <- rep(Q, G)
    if(burnin     < 1)    {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- zinit
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
      z.ind      <- lapply(seq_len(G), function(g) z == g)
      
    # Means
      sum.data   <- do.call(cbind, lapply(seq_len(G), function(g) colSums(data[z.ind[[g]],,drop=F])))
      sum.f      <- lapply(seq_len(G), function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu         <- sim.mu.m(nn=nn, P=P, sigma.mu=sigma.mu, psi.inv=psi.inv,
                             sum.data=sum.data, sum.f=sum.f, lmat=lmat, G=G)
    
    # Scores & Loadings
      c.data     <- lapply(seq_len(G), function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(Q > 0) {
        f        <- sim.score.m(nn=nn, Qs=Qs, lmat=lmat, psi.inv=psi.inv, 
                                c.data=c.data)[obsnames,, drop=F]
        FtF      <- lapply(seq_len(G), function(g) crossprod(f[z.ind[[g]],, drop=F]))
        for(j in seq_len(P)) {
          psi.inv.j <- psi.inv[j,]
          c.data.j  <- lapply(c.data, function(dat) dat[,j])
          lmat[j,,] <- sim.load.m(l.sigma=l.sigma, Qs=Qs, c.data.j=c.data.j, f=f, 
                                  psi.inv.j=psi.inv.j, FtF=FtF, G=G, z.ind=z.ind)
        }
      } else {
        f        <- matrix(, nr=N, nc=0)
        lmat     <- array(, dim=c(P, 0, G))
      }
                  
    # Uniquenesses
      psi.inv    <- sim.psi.im(nn=nn, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                               c.data=c.data, f=f, lmat=lmat, G=G, z.ind=z.ind)
    
    # Mixing Proportions
      pi.prop    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi        <- 1/psi.inv
      Sigma      <- lapply(seq_len(G), function(g) tcrossprod(as.matrix(lmat[,,g])) + diag(psi[,g]))
      z.res      <- sim.z(data=data, mu=mu, Sigma=Sigma, 
                          G=G, P=P, pi.prop=pi.prop)
      z          <- z.res$z
      
      if(is.element(iter, iters)) {
        new.iter <- which(iters == iter)
        log.like <- sum(z.res$log.likes)
        if(sw["mu.sw"])             mu.store[,,new.iter]    <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,new.iter]     <- f
        if(all(sw["l.sw"], Q > 0))  load.store[,,,new.iter] <- lmat
        if(sw["psi.sw"])            psi.store[,,new.iter]   <- psi
        if(sw["pi.sw"])             pi.store[,new.iter]     <- pi.prop
                                    z.store[,new.iter]      <- z 
                                    ll.store[new.iter]      <- log.like
      }  
    }
    returns   <- list(mu       = if(sw["mu.sw"])               mu.store,
                      f        = if(all(sw["f.sw"], Q > 0))    f.store, 
                      load     = if(all(sw["l.sw"], Q > 0))    load.store, 
                      psi      = if(sw["psi.sw"])              psi.store,
                      pi.prop  = if(sw["pi.sw"])               pi.store,
                      z        = z.store,
                      ll.store = ll.store)
    attr(returns, "K")        <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
    attr(returns, "Z.init")   <- zinit
    return(returns)
  }