################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MFA      <- function(Q = NULL, data = NULL, iters = NULL,
                             N = NULL, P = NULL, sigma.mu = NULL,
                             psi.alpha = NULL, psi.beta = NULL, G = NULL,
                             burnin = NULL, thinning = NULL, sw = NULL,
                             verbose = NULL, sigma.l = NULL, alpha.pi = NULL, 
                             zinit = NULL, zlist = NULL,  ...) {
        
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
    pi.alpha     <- rep(alpha.pi, G)
    if(zinit == "priors")  {
      pi.prop    <- sim.pi(pi.alpha=pi.alpha)
      z          <- sim.z.p(N=N, prob.z=pi.prop)
    } else   {
      if(zinit == "list")  {
        z        <- as.numeric(zlist)
      } else {
        z        <- factor(kmeans(data, G, nstart=100)$cluster, levels=Gseq)
      }
      pi.prop    <- prop.table(tabulate(z, nbins=G))
    } 
    zinit        <- z
    if(Q > 0) {
      fact       <- lapply(Gseq, function(g) factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)))
      mu         <- do.call(cbind, lapply(Gseq, function(g) colMeans(data[z == g,, drop=F])))
      f          <- do.call(rbind, lapply(Gseq, function(g) fact[[g]]$scores))
      lmat       <- lapply(Gseq, function(g) fact[[g]]$loadings)
      psi.inv    <- 1/do.call(cbind, lapply(Gseq, function(g) fact[[g]]$uniquenesses))
    } else    {
      mu         <- sim.mu.mp(P=P, mu.sigma=mu.sigma, G=G) 
      f          <- sim.f.mp(N=N, Q=Q)
      lmat       <- sim.load.mp(P=P, Q=Q, l.sigma=l.sigma, G=G)
      psi.inv    <- sim.psi.imp(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, G=G)
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
      z.ind      <- lapply(Gseq, function(g) z == g)
      
    # Means
      sum.data   <- do.call(cbind, lapply(Gseq, function(g) colSums(data[z.ind[[g]],,drop=F])))
      sum.f      <- lapply(Gseq, function(g) colSums(f[z.ind[[g]],, drop=F]))
      mu         <- sim.mu.m(nn=nn, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv,
                             sum.data=sum.data, sum.f=sum.f, lmat=lmat, G=G)
    
    # Scores & Loadings
      c.data     <- lapply(Gseq, function(g) sweep(data[z.ind[[g]],, drop=F], 2, mu[,g], FUN="-"))
      if(Q > 0) {
        f        <- sim.score.m(nn=nn, Qs=Qs, lmat=lmat, psi.inv=psi.inv, 
                                c.data=c.data)[obsnames,, drop=F]
        FtF      <- lapply(Gseq, function(g) crossprod(f[z.ind[[g]],, drop=F]))
        lmat     <- sim.load.m(l.sigma=l.sigma, Qs=Qs, c.data=c.data, f=f, P=P,
                               psi.inv=psi.inv, FtF=FtF, G=G, z.ind=z.ind)
      }
     #} else {
     #  f        <- matrix(, nr=N, nc=0)
     #  lmat     <- lapply(Gseq, function(g) matrix(, nr=P, nc=0))
     #}
                  
    # Uniquenesses
      psi.inv    <- sim.psi.im(nn=nn, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                               c.data=c.data, f=f, lmat=lmat, G=G, z.ind=z.ind)
    
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
        if(sw["mu.sw"])             mu.store[,,new.it]     <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,new.it]      <- f
        if(all(sw["l.sw"], Q > 0)) {
          for(g in Gseq)  {
                                    load.store[,,g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])            psi.store[,,new.it]    <- psi
        if(sw["pi.sw"])             pi.store[,new.it]      <- pi.prop
                                    z.store[,new.it]       <- z 
                                    ll.store[new.it]       <- log.like
      }  
    }
    returns   <- list(mu       = if(sw["mu.sw"])              mu.store,
                      f        = if(all(sw["f.sw"], Q > 0))   f.store, 
                      load     = if(all(sw["l.sw"], Q > 0))   load.store, 
                      psi      = if(sw["psi.sw"])             psi.store,
                      pi.prop  = if(sw["pi.sw"])              pi.store,
                      z        = z.store,
                      ll.store = ll.store)
    attr(returns, "K")        <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
    attr(returns, "Z.init")   <- zinit
    return(returns)
  }