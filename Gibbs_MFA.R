################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MFA        <- function(Q, data, iters, N, P, G, mu.zero,
                               sigma.mu, sigma.l, burnin, mu,
                               thinning, psi.alpha, psi.beta,
                               sw, verbose, cluster, ...) {
         
  # Define & initialise variables
    n.iters        <- round(max(iters), -1)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Gseq)
    iternames      <- paste0("Iteration", seq_len(n.store))
    Q0             <- Q > 0
    Q0s            <- rep(Q0, G)
    Q1             <- Q > 1
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store      <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q0) facnames, gnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store    <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(varnames, gnames, iternames)
    }
    if(sw["pi.sw"])  {
      pi.store     <- matrix(0, nr=G, nc=n.store)
      dimnames(pi.store)   <- list(gnames, iternames)
    }
    z.store        <- matrix(0, nr=N, nc=n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    ll.store       <- rep(0, n.store)
    
    mu.sigma       <- 1/sigma.mu
    if(all(mu.zero == 0)) {
      mu.zero      <- matrix(0, nr=1, nc=G)
      cluster$l.switch[1]  <- F
    }
    l.sigma        <- 1/sigma.l 
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    pi.alpha       <- cluster$pi.alpha
    pi.prop        <- cluster$pi.prop
    mu0g           <- cluster$l.switch[1]
    psi0g          <- cluster$l.switch[2]
    label.switch   <- any(cluster$l.switch)
    f              <- sim.f.p(N=N, Q=Q)
    lmat           <- lapply(Gseq, function(g) sim.load.p(Q=Q, P=P, sigma.l=sigma.l, shrink=F))
    psi.inv        <- do.call(cbind, lapply(Gseq, function(g) sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g])))
    if(Q0) {
      fact.ind     <- nn <= 2.5 * Q
      fail.gs      <- which(fact.ind)
      for(g in which(!fact.ind)) {
        fact       <- try(factanal(data[z == g,, drop=F], factors=Q, scores="regression", control=list(nstart=50)), silent=T)
        if(!inherits(fact, "try-error")) {
          f[z == g,]       <- fact$scores
          lmat[[g]]        <- fact$loadings
          psi.inv[,g]      <- 1/fact$uniquenesses
        } else {
          fail.gs  <- c(fail.gs, g)
        }
      }
      fail.gs      <- fail.gs[order(fail.gs)]
      len.fail     <- length(fail.gs)
      if(len.fail   > 0)      message(paste0("Parameters of the following group", ifelse(len.fail > 2, "s ", " "), "were initialised by simulation from priors, not factanal: ", ifelse(len.fail > 1, paste0(paste0(fail.gs[-len.fail], sep="", collapse=", "), " and ", fail.gs[len.fail]), fail.gs), " - G=", G, ", Q=", Q))
    } else     {
      nn           <- tabulate(z, nbins=G)
      psi.tmp      <- psi.inv
      psi.inv      <- do.call(cbind, lapply(Gseq, function(g) if(nn[g] > 1) 1/apply(data[z == g,, drop=F], 2, var) else psi.tmp[,g]))
      inf.ind      <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    l.sigma        <- l.sigma * diag(Q)
    lmat           <- array(unlist(lmat, use.names=F), dim=c(P, Q, G))
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      f.store[,,1]         <- f
      load.store[,,,1]     <- lmat
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, Sigma=lapply(Gseq,
                                  function(g) tcrossprod(lmat[,,g]) + diag(1/psi.inv[,g])), Q0=Q0s)$log.likes)
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
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      z.ind        <- lapply(Gseq, function(g) z == g)
      dat.g        <- lapply(Gseq, function(g) data[z.ind[[g]],, drop=F])
      
    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(Q0) {
        f.tmp      <- lapply(Gseq, function(g) if(nn0[g]) sim.score(N=nn[g], lmat=lmat[,,g], Q=Q, c.data=c.data[[g]], psi.inv=psi.inv[,g], Q1=Q1) else matrix(, nr=0, nc=Q))
        FtF        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(f.tmp[[g]]))
        lmat       <- array(unlist(lapply(Gseq, function(g) if(nn0[g]) matrix(unlist(lapply(Pseq, function(j) sim.load(l.sigma=l.sigma, Q=Q, c.data=c.data[[g]][,j], 
                            P=P, f=f.tmp[[g]], Q1=Q1, psi.inv=psi.inv[,g][j], FtF=FtF[[g]], shrink=F)), use.names=F), nr=P, byrow=T)), use.names=F), dim=c(P, Q, G))
        f          <- do.call(rbind, f.tmp)[obsnames,, drop=F]
      } else {
        f.tmp      <- lapply(Gseq, function(g) f[z.ind[[g]],, drop=F])
      }
      
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.f        <- lapply(f.tmp, colSums)
      mu           <- do.call(cbind, lapply(Gseq, function(g) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], sum.f=sum.f[[g]],
                              sum.data=sum.data[[g]], P=P, lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g]), mu.zero=mu.zero[,g])))
      
    # Uniquenesses
      psi.inv      <- do.call(cbind, lapply(Gseq, function(g) if(nn0[g]) sim.psi.i(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], P=P, f=f.tmp[[g]],
                              psi.beta=psi.beta[,g], lmat=if(Q1) lmat[,,g] else as.matrix(lmat[,,g])) else sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)))
    
    # Mixing Proportions
      pi.prop      <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi          <- 1/psi.inv
      Sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[,,g]) + diag(psi[,g]))
      z.res        <- sim.z(data=data, mu=mu, Sigma=Sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0s)
      z            <- z.res$z
    
    # Label Switching
      if(label.switch) {
        switch.lab <- lab.switch(z.new=z, z.old=z.temp, Gs=Gseq)
        z          <- switch.lab$z
        z.perm     <- switch.lab$z.perm
        perm       <- identical(as.integer(z.perm), Gseq)
        if(!perm)  {
          mu       <- mu[,z.perm, drop=F]
          lmat     <- lmat[,,z.perm, drop=F]
          psi.inv  <- psi.inv[,z.perm, drop=F]
          pi.prop  <- pi.prop[,z.perm, drop=F]
         if(mu0g)  {
          mu.zero  <- mu.zero[,z.perm, drop=F]
         }
         if(psi0g) {
          psi.beta <- psi.beta[,z.perm, drop=F]
         }
        }
      }
      
      if(is.element(iter, iters))  {
        new.it     <- which(iters == iter)
        log.like   <- sum(z.res$log.likes) 
        if(sw["mu.sw"])            mu.store[,,new.it]      <- mu  
        if(all(sw["f.sw"], Q0))    f.store[,,new.it]       <- f
        if(all(sw["l.sw"], Q0))    load.store[,,,new.it]   <- lmat
        if(sw["psi.sw"])           psi.store[,,new.it]     <- psi
        if(sw["pi.sw"])            pi.store[,new.it]       <- pi.prop
                                   z.store[,new.it]        <- z 
                                   ll.store[new.it]        <- log.like
      }  
    }
    returns        <- list(mu       = if(sw["mu.sw"])         mu.store,
                           f        = if(all(sw["f.sw"], Q0)) f.store, 
                           load     = if(all(sw["l.sw"], Q0)) load.store, 
                           psi      = if(sw["psi.sw"])        psi.store,
                           pi.prop  = if(sw["pi.sw"])         pi.store,
                           z.store  = z.store,
                           ll.store = ll.store)
    attr(returns, "K")  <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
    return(returns)
  }