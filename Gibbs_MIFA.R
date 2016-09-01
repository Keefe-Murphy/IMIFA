################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MIFA       <- function(Q, data, iters, N, P, G, mu.zero,
                               sigma.mu, burnin, thinning, mu,
                               psi.alpha, psi.beta, verbose, 
                               sw, cluster, phi.nu, b0, b1, prop,
                               beta.d1, beta.dk, adapt, epsilon, ...) {
        
  # Define & initialise variables
    start.time     <- proc.time()
    total          <- max(iters)
    if(verbose)       pb   <- txtProgressBar(min=0, max=total, style=3)
    n.store        <- length(iters)
    Gseq           <- seq_len(G)
    Pseq           <- seq_len(P)
    obsnames       <- rownames(data)
    varnames       <- colnames(data)
    colnames(data) <- NULL
    facnames       <- paste0("Factor ", seq_len(Q))
    gnames         <- paste0("Group ", Gseq)
    iternames      <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store     <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(varnames, gnames, iternames)
    }
    if(sw["s.sw"])   {
      eta.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)  <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store   <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, gnames, iternames)
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
    ll.store       <- rep(0, n.store)
    Q.star         <- Q
    Qs             <- rep(Q, G)
    Q.store        <- matrix(0, nr=G, nc=n.store)
    dimnames(z.store)      <- list(obsnames, iternames)
    dimnames(Q.store)      <- list(gnames, iternames)
    
    mu.sigma       <- 1/sigma.mu
    if(all(mu.zero == 0)) {
      mu.zero      <- matrix(0, nr=1, nc=G)
      cluster$l.switch[1]  <- FALSE
    }
    z              <- cluster$z
    z.temp         <- factor(z, levels=Gseq)
    nn             <- tabulate(z, nbins=G)
    pi.alpha       <- cluster$pi.alpha
    pi.prop        <- cluster$pi.prop
    alpha.d1       <- cluster$alpha.d1
    alpha.dk       <- cluster$alpha.dk
    ad1.x          <- length(unique(alpha.d1)) == 1
    adk.x          <- length(unique(alpha.dk)) == 1
    mu0g           <- cluster$l.switch[1]
    psi0g          <- cluster$l.switch[2]
    delta0g        <- cluster$l.switch[3]
    qstar0g        <- cluster$l.switch[4]
    label.switch   <- any(cluster$l.switch)
    eta            <- sim.eta.p(N=N, Q=Q)
    phi            <- lapply(Gseq, function(g) sim.phi.p(Q=Q, P=P, phi.nu=phi.nu))
    delta          <- lapply(Gseq, function(g) c(sim.delta.p(alpha=alpha.d1[g], beta=beta.d1), sim.delta.p(Q=Q, alpha=alpha.dk[g], beta=beta.dk)))
    tau            <- lapply(delta, cumprod)
    lmat           <- lapply(Gseq, function(g) matrix(unlist(lapply(Pseq, function(j) sim.load.p(Q=Q, phi=phi[[g]][j,], tau=tau[[g]], P=P)), use.names=FALSE), nr=P, byrow=TRUE))
    psi.inv        <- vapply(Gseq, function(g) sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    fact.ind       <- nn <= 2.5 * Q
    fail.gs        <- which(fact.ind)
    for(g in which(!fact.ind)) {
      fact         <- try(factanal(data[z == g,, drop=FALSE], factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
      if(!inherits(fact, "try-error")) {
        eta[z == g,]       <- fact$scores
        lmat[[g]]          <- fact$loadings
        psi.inv[,g]        <- 1/fact$uniquenesses
      } else  {
        fail.gs    <- c(fail.gs, g)
      }               
    }
    fail.gs        <- fail.gs[order(fail.gs)]
    len.fail       <- length(fail.gs)
    if(len.fail     > 0)      message(paste0("Parameters of the following group", ifelse(len.fail > 2, "s ", " "), "were initialised by simulation from priors, not factanal: ", ifelse(len.fail > 1, paste0(paste0(fail.gs[-len.fail], sep="", collapse=", "), " and ", fail.gs[len.fail]), fail.gs), " - G=", G, ", Q=", Q))
    if(burnin       < 1)  {
      mu.store[,,1]        <- mu
      eta.store[,,1]       <- eta
      load.store[,,,1]     <- array(unlist(lmat, use.names=FALSE), dim=c(P, Q, G))
      psi.store[,,1]       <- 1/psi.inv
      pi.store[,1]         <- pi.prop
      z.store[,1]          <- z
      ll.store[1]          <- sum(sim.z(data=data, mu=mu, Gseq=Gseq, N=N, pi.prop=pi.prop, sigma=lapply(Gseq,
                                  function(g) tcrossprod(lmat[[g]]) + diag(1/psi.inv[,g])), Q0=Qs > 0)$log.likes)
      Q.store[,1]          <- Qs
    }
    init.time      <- proc.time() - start.time
    
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose   && iter  < burnin) setTxtProgressBar(pb, iter)
      
    # Mixing Proportions
      pi.prop[]    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
      
    # Cluster Labels
      psi          <- 1/psi.inv
      sigma        <- lapply(Gseq, function(g) tcrossprod(lmat[[g]]) + diag(psi[,g]))
      Q0           <- Qs > 0
      Q1           <- Qs > 1
      z.res        <- sim.z(data=data, mu=mu, sigma=sigma, Gseq=Gseq, N=N, pi.prop=pi.prop, Q0=Q0)
      z            <- z.res$z
      nn           <- tabulate(z, nbins=G)
      nn0          <- nn > 0
      z.ind        <- lapply(Gseq, function(g) z == g)
      dat.g        <- lapply(Gseq, function(g) data[z.ind[[g]],, drop=FALSE])
    
    # Scores & Loadings
      c.data       <- lapply(Gseq, function(g) sweep(dat.g[[g]], 2, mu[,g], FUN="-"))
      if(!any(Q0))    {
        eta        <- matrix(, nr=N, nc=0)
        eta.tmp    <- lapply(Gseq, function(g) eta[z.ind[[g]],, drop=FALSE])
        lmat       <- lapply(Gseq, function(g) matrix(, nr=P, nc=0))
      } else {
        eta.tmp    <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) sim.score(N=nn[g], lmat=lmat[[g]], Q=Qs[g], Q1=Q1[g], c.data=c.data[[g]], psi.inv=psi.inv[,g]) else matrix(, nr=ifelse(Q0[g], 0, nn[g]), nc=Qs[g]))
        EtE        <- lapply(Gseq, function(g) if(nn0[g]) crossprod(eta.tmp[[g]]))
        lmat       <- lapply(Gseq, function(g) if(all(nn0[g], Q0[g])) matrix(unlist(lapply(Pseq, function(j) sim.load(Q=Qs[g], P=P, c.data=c.data[[g]][,j], EtE=EtE[[g]],
                             eta=eta.tmp[[g]], psi.inv=psi.inv[,g][j], Q1=Q1[g], phi=phi[[g]][j,], tau=tau[[g]], shrink=TRUE)), use.names=FALSE), nr=P, byrow=TRUE) else 
                             matrix(unlist(lapply(Pseq, function(j) sim.load.p(Q=Qs[g], phi=phi[[g]][j,], tau=tau[[g]], P=P)), use.names=FALSE), nr=P, byrow=FALSE))
        eta.tmp    <- if(length(unique(Qs)) != 1) lapply(Gseq, function(g) cbind(eta.tmp[[g]], matrix(0, nr=nn[g], nc=max(Qs) - Qs[g]))) else eta.tmp
        q0ng       <- !Q0 & nn0
        if(any(q0ng)) {
          eta.tmp[q0ng]    <- lapply(Gseq[q0ng], function(g, x=eta.tmp[[g]]) { row.names(x) <- obsnames[z.ind[[g]]]; x })
        }
        eta        <- do.call(rbind, eta.tmp)[obsnames,, drop=FALSE]
      }
      
    # Means
      sum.data     <- lapply(dat.g, colSums)
      sum.eta      <- lapply(eta.tmp, colSums)
      mu           <- vapply(Gseq, function(g) if(nn0[g]) sim.mu(N=nn[g], mu.sigma=mu.sigma, psi.inv=psi.inv[,g], P=P, sum.eta=sum.eta[[g]][seq_len(Qs[g])], 
                             sum.data=sum.data[[g]], lmat=lmat[[g]], mu.zero=mu.zero[,g]) else sim.mu.p(P=P, sigma.mu=sigma.mu, mu.zero=mu.zero[,g]), numeric(P))
      
    # Uniquenesses
      psi.inv      <- vapply(Gseq, function(g) if(nn0[g]) sim.psi.inv(N=nn[g], psi.alpha=psi.alpha, c.data=c.data[[g]], psi.beta=psi.beta[,g], lmat=lmat[[g]],
                             P=P, eta=eta.tmp[[g]][,seq_len(Qs[g]), drop=FALSE]) else sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta[,g]), numeric(P))
    
    # Local Shrinkage
      load.2       <- lapply(lmat, function(lg) lg * lg)
      phi          <- lapply(Gseq, function(g) sim.phi(Q=Qs[g], P=P, phi.nu=phi.nu, tau=tau[[g]], load.2=load.2[[g]]))
    
    # Global Shrinkage
      sum.terms    <- lapply(Gseq, function(g) diag(crossprod(phi[[g]], load.2[[g]])))
      for(g in Gseq)  {
        Qg         <- Qs[g]
        for(k in seq_len(Qg)) { 
          delta[[g]][k]    <- if(k > 1) sim.deltak(alpha.dk=alpha.dk[g], beta.dk=beta.dk, delta.k=delta[[g]][k], Q=Qg, P=P, 
                              k=k, tau.kq=tau[[g]][k:Qg], sum.term.kq=sum.terms[[g]][k:Qg]) else sim.delta1(alpha.d1=alpha.d1[g],
                              beta.d1=beta.d1, delta.1=delta[[g]][1], Q=Qg, P=P, tau=tau[[g]], sum.term=sum.terms[[g]])
          tau[[g]]         <- cumprod(delta[[g]])
        }
      }
    
    # Adaptation  
      if(all(adapt, iter > burnin)) {      
        if(runif(1) < 1/exp(b0 + b1 * (iter - burnin))) {
          lind     <- lapply(Gseq, function(g) if(Q0[g]) colSums(abs(lmat[[g]]) < epsilon)/P else 0)
          colvec   <- lapply(lind, function(lx) lx >= prop)
          nonred   <- lapply(colvec, function(cv) which(cv == 0))
          numred   <- lapply(colvec, sum)
          notred   <- vapply(Gseq, function(g) numred[[g]] == 0, logical(1))
          Qs.old   <- Qs
          Qs       <- vapply(Gseq, function(g) if(notred[g]) Qs.old[g] + 1 else Qs.old[g] - numred[[g]], numeric(1))
          phi      <- lapply(Gseq, function(g) if(notred[g]) cbind(phi[[g]][,seq_len(Qs.old[g])], rgamma(n=P, shape=phi.nu, rate=phi.nu)) else phi[[g]][,nonred[[g]], drop=FALSE])
          delta    <- lapply(Gseq, function(g) if(notred[g]) c(delta[[g]][seq_len(Qs.old[g])], rgamma(n=1, shape=alpha.dk[g], rate=beta.dk)) else delta[[g]][nonred[[g]]])  
          tau      <- lapply(delta, cumprod)
          lmat     <- lapply(Gseq, function(g, Qg=Qs[g]) if(notred[g]) cbind(lmat[[g]][,seq_len(Qs.old[g])], rnorm(n=P, mean=0, sd=sqrt(1/(phi[[g]][,Qg] * tau[[g]][Qg])))) else lmat[[g]][,nonred[[g]], drop=FALSE])
          eta      <- if(max(Qs) > max(Qs.old)) cbind(eta[,seq_len(max(Qs.old))], rnorm(N)) else eta[,seq_len(max(Qs)), drop=FALSE]
        }
      }
    
    # Label Switching
      if(label.switch)   {
        switch.lab <- lab.switch(z.new=z, z.old=z.temp, Gs=Gseq)
        z          <- switch.lab$z
        z.perm     <- switch.lab$z.perm
        if(!identical(as.integer(z.perm), Gseq)) {
         if(length(unique(Qs)) != 1) {
          Qs       <- Qs[z.perm]  
         }
         mu        <- mu[,z.perm, drop=FALSE]
         lmat      <- lmat[z.perm]
         delta     <- delta[z.perm]
         phi       <- phi[z.perm]
         tau       <- tau[z.perm]
         psi.inv   <- psi.inv[,z.perm, drop=FALSE]
         pi.prop   <- pi.prop[z.perm]
         if(mu0g)        {
          mu.zero  <- mu.zero[,z.perm, drop=FALSE]
         }
         if(psi0g)       {
          psi.beta <- psi.beta[,z.perm, drop=FALSE]
         }
         if(all(delta0g, 
                !ad1.x)) {
          alpha.d1 <- alpha.d1[z.perm]
         }
         if(all(delta0g, 
                !adk.x)) {
          alpha.dk <- alpha.dk[z.perm]
         } 
        }
      }
      
    if(any(Qs > Q.star))      stop(paste0("Q cannot exceed initial number of loadings columns: try increasing range.Q from ", Q.star))
      if(is.element(iter, iters))    {
        if(verbose)   setTxtProgressBar(pb, iter)
        new.it     <- which(iters == iter)
        if(sw["mu.sw"])    mu.store[,,new.it]       <- mu  
        if(all(sw["s.sw"], 
           any(Q0)))  eta.store[,seq_len(max(Qs)),new.it]  <- eta
        if(sw["l.sw"])   {
          for(g in Gseq) {
            if(Q0[g]) load.store[,seq_len(Qs[g]),g,new.it] <- lmat[[g]]
          }
        }
        if(sw["psi.sw"])   psi.store[,,new.it]      <- psi
        if(sw["pi.sw"])    pi.store[,new.it]        <- pi.prop
                           z.store[,new.it]         <- z 
                           ll.store[new.it]         <- sum(z.res$log.likes)
                           Q.store[,new.it]         <- Qs
      }
    }
    close(pb)
    Qmax           <- seq_len(max(Q.store))
    returns        <- list(mu       = if(sw["mu.sw"])  mu.store,
                           eta      = if(sw["s.sw"])   as.simple_sparse_array(eta.store[,Qmax,, drop=FALSE]), 
                           load     = if(sw["l.sw"])   as.simple_sparse_array(load.store[,Qmax,,, drop=FALSE]), 
                           psi      = if(sw["psi.sw"]) psi.store,
                           pi.prop  = if(sw["pi.sw"])  pi.store,
                           z.store  = z.store,
                           ll.store = ll.store,
                           Q.store  = Q.store,
                           time     = init.time)
    return(returns)
  }