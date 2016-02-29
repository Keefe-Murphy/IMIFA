################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Gibbs Sampler Function
  gibbs.FA       <- function(Q, data, n.iters, N, P, 
                             sigma.mu, psi.alpha, psi.beta,
                             burnin, thinning, n.store,
                             verbose, sw, sigma.l, ...) {
        
  # Define & initialise variables
    cnames       <- colnames(data)
    n.store      <- ceiling((n.iters - burnin)/thinning)
    if(sw["mu.sw"]) {
      mu.store   <- matrix(0, nr=P, nc=n.store)
      dimnames(mu.store)[[1]]   <- cnames
      dimnames(mu.store)[[2]]   <- paste0("Iteration", 1:n.store)
    }
    if(sw["f.sw"])  {
      f.store    <- array(0, dim=c(N, Q, n.store))
      dimnames(f.store)[[1]]    <- rownames(data)
      if(Q > 0) dimnames(f.store)[[2]]    <- paste0("Factor ", 1:Q)
      dimnames(f.store)[[3]]    <- paste0("Iteration", 1:n.store)
    }
    if(sw["l.sw"])  {
      load.store <- array(0, dim=c(P, Q, n.store))
      dimnames(load.store)[[1]] <- cnames
      if(Q > 0) dimnames(load.store)[[2]] <- paste0("Factor ", 1:Q)
      dimnames(load.store)[[3]] <- paste0("Iteration", 1:n.store)
    }
    if(sw["p.sw"])  {
      psi.store  <- matrix(0, nr=P, nc=n.store)
      dimnames(psi.store)[[1]]  <- cnames
      dimnames(psi.store)[[2]]  <- paste0("Iteration", 1:n.store)
    }
    post.Sigma   <- matrix(0, nr=P, nc=P)
    cov.emp      <- cov(data)
    dimnames(post.Sigma)        <- list(cnames, cnames)
    dimnames(cov.emp)           <- dimnames(post.Sigma)
    
    sigma.mu     <- 1/sigma.mu
    sigma.l      <- 1/sigma.l
    mu           <- sim.mu.p(P, sigma.mu)  
    f            <- sim.f.p(Q, N)
    lmat         <- sim.l.p(Q, P, sigma.l)
    psi.inv      <- sim.pi.p(P, psi.alpha, psi.beta)
    l.sigma      <- sigma.l * diag(Q)
    sum.data     <- colSums(data)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(verbose) {
        if(iter < burnin && iter %% (burnin/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if (iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
    # Means
      sum.f      <- colSums(f)
      mu         <- sim.mu(N, P, sigma.mu, psi.inv, sum.data, sum.f, lmat)
    
    # Scores
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q > 0) {
        f        <- sim.scores(N, Q, lmat, psi.inv, c.data)
      } else {
        f        <- matrix(, nr=N, nc=0)
      }
                
    # Loadings
      FtF        <- crossprod(f)
      if(Q > 0) {
        for (j in 1:P) {
          psi.inv.j <- psi.inv[j]
          c.data.j  <- c.data[,j]
          lmat[j,]  <- sim.load(l.sigma, Q, c.data.j, f, psi.inv.j, FtF)
        }
      } else {
        lmat     <- matrix(, nr=P, nc=0)
      }
      
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N, P, psi.alpha, psi.beta, c.data, f, lmat)
    
      if(iter >= burnin && iter %% thinning == 0) {
        new.iter <- ceiling((iter - burnin)/thinning)
        psi      <- 1/psi.inv
        if(sw["mu.sw"]) mu.store[,new.iter]    <- mu  
        if(sw["f.sw"])  f.store[,,new.iter]    <- f
        if(sw["l.sw"])  load.store[,,new.iter] <- lmat
        if(sw["p.sw"])  psi.store[,new.iter]   <- psi
        Sigma       <-  tcrossprod(lmat) + diag(psi)
        post.Sigma  <-  post.Sigma + Sigma/n.store
      }  
    }
    returns   <- list(mu   = if(sw["mu.sw"]) mu.store,
                      f    = if(sw["f.sw"])  f.store, 
                      load = if(sw["l.sw"])  load.store, 
                      psi  = if(sw["p.sw"])  psi.store,
                      cov.mat    = cov.emp,
                      post.Sigma = post.Sigma)
    return(returns[!sapply(returns, is.null)])
  }