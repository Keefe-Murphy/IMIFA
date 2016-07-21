###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Gibbs Sampler Function
  gibbs.IFA      <- function(Q, data, iters, N, P, sigma.mu, 
                             psi.alpha, psi.beta, burnin, mu,
                             thinning, verbose, sw, mu.zero,
                             phi.nu, alpha.d1, alpha.dk, 
                             beta.d1, beta.dk, b0, b1,
                             adapt, prop, epsilon, ...) {    
    
  # Define & initialise variables
    n.iters      <- round(max(iters), -1)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    iternames    <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"])  {
      mu.store   <- matrix(0, nr=P, nc=n.store)
      dimnames(mu.store)   <- list(varnames, iternames)
    }
    if(sw["f.sw"])   {
      f.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)    <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store  <- matrix(0, nr=P, nc=n.store)
      dimnames(psi.store)  <- list(varnames, iternames)
    }
    post.mu      <- setNames(rep(0, P), varnames)
    post.psi     <- setNames(rep(0, P), varnames)
    cov.emp      <- cov(data)
    cov.est      <- matrix(0, nr=P, nc=P)
    ll.store     <- rep(0, n.store)
    Q.star       <- Q
    Q.store      <- setNames(rep(0, n.store), iternames)
    dimnames(cov.emp)      <- list(varnames, varnames)
    dimnames(cov.est)      <- dimnames(cov.emp)
    
    mu.sigma     <- 1/sigma.mu
    f            <- sim.f.p(Q=Q, N=N)
    psi.inv      <- sim.psi.ip(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    phi          <- sim.phi.p(Q=Q, P=P, phi.nu=phi.nu)
    delta        <- c(sim.delta.p(alpha=alpha.d1, beta=beta.d1), sim.delta.p(Q=Q, alpha=alpha.dk, beta=beta.dk))
    tau          <- cumprod(delta)
    lmat         <- matrix(unlist(lapply(Pseq, function(j) sim.load.p(Q=Q, phi=phi[j,], tau=tau, P=P)), use.names=F), nr=P, byrow=T)
    sum.data     <- mu * N
    if(burnin     < 1) {
      mu.store[,1]         <- mu
      f.store[,,1]         <- f
      load.store[,,1]      <- lmat
      psi.store[,1]        <- 1/psi.inv
      ll.store[1]          <- sum(dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + diag(1/psi.inv), log=T))
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
      Q0         <- Q > 0
      Q1         <- Q > 1
      
    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q0) {
        f        <- sim.score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        FtF      <- crossprod(f)
        lmat     <- matrix(unlist(lapply(Pseq, function(j) sim.load(Q=Q, tau=tau, P=P, f=f, Q1=Q1,
                           c.data=c.data[,j], phi=phi[j,], psi.inv=psi.inv[j], FtF=FtF)), use.names=F), nr=P, byrow=T)
      } else {
        f        <- matrix(, nr=N, nc=0)
        lmat     <- matrix(, nr=P, nc=0)
      }     
      
    # Means
      sum.f      <- colSums(f)
      mu         <- sim.mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.f=sum.f, lmat=lmat, mu.zero=mu.zero)
    
    # Uniquenesses
      psi.inv    <- sim.psi.i(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, c.data=c.data, f=f, lmat=lmat)
    
    # Local Shrinkage
      load.2     <- lmat * lmat
      phi        <- sim.phi(Q=Q, P=P, phi.nu=phi.nu, tau=tau, load.2=load.2)
          
    # Global Shrinkage
      sum.term   <- diag(crossprod(phi, load.2))
    if(Q0) {
      delta[1]   <- sim.delta1(Q=Q, P=P, alpha.d1=alpha.d1, delta=delta,
                               beta.d1=beta.d1, tau=tau, sum.term=sum.term)
      tau        <- cumprod(delta)  
    } 
    if(Q1) {
      for(k in seq_len(Q)[-1]) { 
        delta[k] <- sim.deltak(Q=Q, P=P, k=k, alpha.dk=alpha.dk, delta=delta,
                               beta.dk=beta.dk, tau=tau, sum.term=sum.term)
        tau      <- cumprod(delta)      
      }
    }
    
    # Adaptation  
      if(all(adapt, iter > burnin)) {      
        prob     <- 1/exp(b0 + b1 * (iter - burnin))
        unif     <- runif(n=1, min=0, max=1)     
        if(unif   < prob) { # check whether to adapt or not
          if(Q0) {
            lind <- colSums(abs(lmat) < epsilon) / P
          } else {
            lind <- 0
          }
          colvec <- lind >= prop
          numred <- sum(colvec)
          if(numred == 0) { # simulate extra columns from priors
            Q       <- Q + 1
            f       <- cbind(f, rnorm(n=N, mean=0, sd=1))         
            phi     <- cbind(phi, rgamma(n=P, shape=phi.nu, rate=phi.nu))
            delta   <- c(delta, rgamma(n=1, shape=alpha.dk, rate=beta.dk))
            tau     <- cumprod(delta)
            lmat    <- cbind(lmat, rnorm(n=P, mean=0, sd=sqrt(1/(phi[,Q] * tau[Q]))))
          } else          { # remove redundant columns
            nonred  <- which(colvec == 0)
            Q       <- Q - numred
            f       <- f[,nonred, drop=F]
            phi     <- phi[,nonred, drop=F]
            delta   <- delta[nonred]
            tau     <- cumprod(delta)
            lmat    <- lmat[,nonred, drop=F]
          }
        }
      } 
    
    if(Q > Q.star)  stop(paste0("Q cannot exceed initial number of loadings columns: try increasing range.Q from ", Q.star))
     if(is.element(iter, iters)) {
        new.it   <- which(iters == iter)  
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        Sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + Sigma/n.store
        log.like <- sum(dmvn(X=data, mu=mu, sigma=Sigma, log=T))
        if(sw["mu.sw"])             mu.store[,new.it]              <- mu  
        if(all(sw["f.sw"], Q0))     f.store[,seq_len(Q),new.it]    <- f
        if(all(sw["l.sw"], Q0))     load.store[,seq_len(Q),new.it] <- lmat
        if(sw["psi.sw"])            psi.store[,new.it]             <- psi
                                    Q.store[new.it]                <- Q
                                    ll.store[new.it]               <- log.like
      }
    }
    returns      <- list(mu       = if(sw["mu.sw"])  mu.store,
                         f        = if(sw["f.sw"])   as.simple_sparse_array(f.store), 
                         load     = if(sw["l.sw"])   as.simple_sparse_array(load.store), 
                         psi      = if(sw["psi.sw"]) psi.store,
                         post.mu  = post.mu,
                         post.psi = post.psi,
                         cov.emp  = cov.emp,
                         cov.est  = cov.est,
                         ll.store = ll.store,
                         Q.store  = matrix(Q.store, nr=1))
    return(returns)
  }