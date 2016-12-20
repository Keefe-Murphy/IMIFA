###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Gibbs Sampler Function
  .gibbs_IFA     <- function(Q, data, iters, N, P, sigma.mu, sw, mu, prop, 
                             psi.alpha, psi.beta, burnin, thinning, verbose,
                             epsilon, mu.zero, nu, adapt, adapt.at, b0, b1,
                             alpha.d1, alpha.d2, beta.d1, beta.d2, nuplus1, ...) {    
    
  # Define & initialise variables
    start.time   <- proc.time()
    total        <- max(iters)
    if(verbose)     pb     <- txtProgressBar(min=0, max=total, style=3)
    n.store      <- length(iters)
    Pseq         <- seq_len(P)
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    iternames    <- paste0("Iteration", seq_len(n.store))
    dimnames(data)         <- NULL
    if(sw["mu.sw"])  {
      mu.store   <- provideDimnames(matrix(0, nr=P, nc=n.store), base=list(varnames, iternames))
    }
    if(sw["s.sw"])   {
      eta.store  <- array(0, dim=c(N, Q, n.store))
      dimnames(eta.store)  <- list(obsnames, if(Q > 0) facnames, iternames)
    }
    if(sw["l.sw"])   {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store) <- list(varnames, if(Q > 0) facnames, iternames)
    }
    if(sw["psi.sw"]) {
      psi.store  <- provideDimnames(matrix(0, nr=P, nc=n.store), base=list(varnames, iternames))
    }
    post.mu      <- setNames(rep(0, P), varnames)
    post.psi     <- setNames(rep(0, P), varnames)
    cov.emp      <- provideDimnames(cova(as.matrix(data)), base=list(varnames, varnames))
    cov.est      <- provideDimnames(matrix(0, nr=P, nc=P), base=dimnames(cov.emp))
    ll.store     <- setNames(rep(0, n.store), iternames)
    Q.star       <- Q
    Q.store      <- setNames(rep(0, n.store), iternames)
    Q.large      <- Q.big  <- FALSE
    
    mu.sigma     <- 1/sigma.mu
    psi.beta     <- unique(round(psi.beta, min(nchar(psi.beta))))
    eta          <- .sim.eta.p(Q=Q, N=N)
    phi          <- .sim.phi.p(Q=Q, P=P, nu=nu, plus1=nuplus1)
    delta        <- c(.sim.delta.p(alpha=alpha.d1, beta=beta.d1), .sim.delta.p(Q=Q, alpha=alpha.d2, beta=beta.d2))
    tau          <- cumprod(delta)
    lmat         <- matrix(unlist(lapply(Pseq, function(j) .sim.load.ps(Q=Q, phi=phi[j,], tau=tau)), use.names=FALSE), nr=P, byrow=TRUE)
    psi.inv      <- .sim.psi.i.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    if(all(Q  < P - sqrt(P + Q), N > P)) {
      fact       <- try(factanal(data, factors=Q, scores="regression", control=list(nstart=50)), silent=TRUE)
      if(!inherits(fact, "try-error"))   {
        eta      <- fact$scores
        lmat     <- fact$loadings
        psi.inv  <- 1/fact$uniquenesses
      }
    } else {
      psi.tmp    <- psi.inv
      psi.inv    <- 1/colVars(data)
      inf.ind    <- is.infinite(psi.inv)
      psi.inv[inf.ind]     <- psi.tmp[inf.ind]
    }
    sum.data     <- mu * N
    if(burnin     < 1) {
      mu.store[,1]         <- mu
      eta.store[,,1]       <- eta
      load.store[,,1]      <- lmat
      psi.store[,1]        <- 1/psi.inv
      ll.store[1]          <- sum(dmvn(X=data, mu=mu, sigma=tcrossprod(lmat) + diag(1/psi.inv), log=TRUE))
    }
    init.time    <- proc.time() - start.time
  
  # Iterate
    for(iter in seq_len(total)[-1]) { 
      if(verbose && iter    < burnin) setTxtProgressBar(pb, iter)
      Q0         <- Q  > 0
      Q1         <- Q == 1
      
    # Scores & Loadings
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q0) {
        eta      <- .sim.score(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data, Q1=Q1)
        lmat     <- matrix(unlist(lapply(Pseq, function(j) .sim.load.s(Q=Q, tau=tau, eta=eta, c.data=c.data[,j], Q1=Q1, 
                           phi=phi[j,], psi.inv=psi.inv[j], EtE=crossprod(eta))), use.names=FALSE), nr=P, byrow=TRUE)
      } else {
        eta      <- base::matrix(0, nr=N, nc=0)
        lmat     <- base::matrix(0, nr=P, nc=0)
      }     
      
    # Means
      mu[]       <- .sim.mu(N=N, P=P, mu.sigma=mu.sigma, psi.inv=psi.inv, sum.data=sum.data, sum.eta=colSums(eta), lmat=lmat, mu.zero=mu.zero)
    
    # Uniquenesses
      psi.inv    <- .sim.psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, c.data=c.data, eta=eta, lmat=lmat)
    
    # Local Shrinkage
      load.2     <- lmat * lmat
      phi        <- .sim.phi(Q=Q, P=P, nu=nu, tau=tau, load.2=load.2, plus1=nuplus1)
          
    # Global Shrinkage
      sum.term   <- diag(crossprod(phi, load.2))
      for(k in seq_len(Q)) { 
        delta[k] <- if(k > 1) .sim.deltak(alpha.d2=alpha.d2, beta.d2=beta.d2, delta.k=delta[k], Q=Q, P=P, k=k,
                    tau.kq=tau[k:Q], sum.term.kq=sum.term[k:Q]) else .sim.delta1(Q=Q, P=P, tau=tau, sum.term=sum.term,
                    alpha.d1=ifelse(Q1, alpha.d2, alpha.d1), beta.d1=ifelse(Q1, beta.d2, beta.d1), delta.1=delta[1])
        tau      <- cumprod(delta)      
      }
    
    # Adaptation  
      if(all(adapt, iter > adapt.at)) {      
        if(runif(1)  < ifelse(iter < burnin, 0.5, 1/exp(b0 + b1 * (iter - adapt.at)))) {
          colvec <- (if(Q0) colSums(abs(lmat) < epsilon) / P else 0) >= prop
          numred <- sum(colvec)
          if(numred == 0) { # simulate extra columns from priors
            Q    <- Q + 1
            Q.big   <- Q > Q.star
            if(Q.big) {
              Q     <- Q.star
            } else {
              eta   <- cbind(eta, rnorm(N))         
              phi   <- cbind(phi, rgamma(n=P, shape=nu + nuplus1, rate=nu))
              delta <- c(delta, rgamma(n=1, shape=alpha.d2, rate=beta.d2))
              tau   <- cumprod(delta)
              lmat  <- cbind(lmat, rnorm(n=P, mean=0, sd=sqrt(1/(phi[,Q] * tau[Q]))))  
            }
          } else          { # remove redundant columns
            nonred  <- which(colvec == 0)
            Q       <- Q - numred
            eta     <- eta[,nonred, drop=FALSE]
            phi     <- phi[,nonred, drop=FALSE]
            delta   <- delta[nonred]
            tau     <- cumprod(delta)
            lmat    <- lmat[,nonred, drop=FALSE]
          }
        }
      } 
    
      if(Q.big && !Q.large && iter > burnin) {       warning(paste0("Q has exceeded initial number of loadings columns since burnin: consider increasing range.Q from ", Q.star), call.=FALSE)
        Q.large  <- TRUE
      }
      if(is.element(iter, iters))  {
        if(verbose) setTxtProgressBar(pb, iter)
        new.it   <- which(iters == iter)  
        psi      <- 1/psi.inv
        post.mu  <- post.mu + mu/n.store
        post.psi <- post.psi + psi/n.store
        sigma    <- tcrossprod(lmat) + diag(psi)
        cov.est  <- cov.est + sigma/n.store
        if(sw["mu.sw"])             mu.store[,new.it]              <- mu  
        if(all(sw["s.sw"], Q0))     eta.store[,seq_len(Q),new.it]  <- eta
        if(all(sw["l.sw"], Q0))     load.store[,seq_len(Q),new.it] <- lmat
        if(sw["psi.sw"])            psi.store[,new.it]             <- psi
                                    Q.store[new.it]                <- Q
                                    ll.store[new.it]               <- sum(dmvn(X=data, mu=mu, sigma=sigma, log=TRUE))
      }
    }
    close(pb)
    Qmax         <- seq_len(max(Q.store))
    returns      <- list(mu       = if(sw["mu.sw"])  mu.store,
                         eta      = if(sw["s.sw"])   as.simple_sparse_array(eta.store[,Qmax,, drop=FALSE]), 
                         load     = if(sw["l.sw"])   as.simple_sparse_array(load.store[,Qmax,, drop=FALSE]), 
                         psi      = if(sw["psi.sw"]) psi.store,
                         post.mu  = post.mu,
                         post.psi = post.psi,
                         cov.emp  = cov.emp,
                         cov.est  = cov.est,
                         ll.store = ll.store,
                         Q.store  = matrix(Q.store, nr=1),
                         time     = init.time)
    attr(returns, "Q.big") <- Q.large
    return(returns)
  }