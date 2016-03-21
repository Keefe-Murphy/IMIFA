################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.FA       <- function(Q = NULL, data = NULL, n.iters = NULL,
                             N = NULL, P = NULL, sigma.mu = NULL,
                             psi.alpha = NULL, psi.beta = NULL,
                             burnin = NULL, thinning = NULL, sw = NULL,
                             n.store = NULL, verbose = NULL, sigma.l = NULL,
                             alpha.pi = NULL, z.init = NULL, z.list = NULL,  ...) {
        
  # Define & initialise variables
    cnames       <- colnames(data)
    rnames       <- rownames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    gnames       <- paste0("Group ", seq_len(G))
    iternames    <- paste0("Iteration", seq_len(n.store))
    if(sw["mu.sw"]) {
      mu.store   <- array(0, dim=c(P, G, n.store))
      dimnames(mu.store)   <- list(cnames, gnames, iternames)
    }
    if(sw["f.sw"])  {
      f.store    <- array(0, dim=c(N, Q, G, n.store))
      dimnames(f.store)    <- list(rnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["l.sw"])  {
      load.store <- array(0, dim=c(P, Q, G, n.store))
      dimnames(load.store) <- list(cnames, if(Q > 0) facnames, gnames, iternames)
    }
    if(sw["p.sw"])  {
      psi.store  <- array(0, dim=c(P, G, n.store))
      dimnames(psi.store)  <- list(cnames, gnames, iternames)
    }
    post.mu      <- matrix(0, nr=P, nc=G)
    post.psi     <- matrix(0, nr=P, nc=G)
   #post.Sigma   <- array(0, dim=c(P, P, G))
    bic.mcmc     <- - Inf
    K            <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
    pen          <- K * log(N)
   #cov.emp      <- cov(data)
    dimnames(post.mu)      <- dimnames(post.psi)   <- list(cnames, gnames)
   #dimnames(post.Sigma)   <- list(cnames, cnames, gnames)
   #dimnames(cov.emp)      <- dimnames(post.Sigma)
    
    sigma.mu     <- 1/sigma.mu
    sigma.l      <- 1/sigma.l
    mu           <- sim.mu.p(P=P, sigma.mu=sigma.mu)  
    f            <- sim.f.p(Q=Q, N=N)
    lmat         <- sim.l.p(Q=Q, P=P, sigma.l=sigma.l)
    psi.inv      <- sim.pi.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta)
    alpha.pi     <- rep(alpha.pi, G)
    pi.prop      <- sim.pi(alpha.pi=alpha.pi)
    if(z.init == "list") {
      z          <- z.list[[G]]
    } else if(z.init == "kmeans") {
      z          <- kmeans(data, G, nstart=100)$cluster
    } else {
      z          <- sim.z.p(N=N, prob.z=pi.prop)
    }
    l.sigma      <- sigma.l * diag(Q)
    sum.data     <- colSums(data)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(verbose) {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
    # Means
      sum.f      <- colSums(f)
      mu         <- sim.mu(N=N, P=P, sigma.mu=sigma.mu, psi.inv=psi.inv,
                           sum.data=sum.data, sum.f=sum.f, lmat=lmat)
    
    # Scores
      c.data     <- sweep(data, 2, mu, FUN="-")
      if(Q > 0) {
        f        <- sim.scores(N=N, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data)
      } else {
        f        <- matrix(, nr=N, nc=0)
      }
                
    # Loadings
      FtF        <- crossprod(f)
      if(Q > 0) {
        for (j in 1:P) {
          psi.inv.j <- psi.inv[j]
          c.data.j  <- c.data[,j]
          lmat[j,]  <- sim.load(l.sigma=l.sigma, Q=Q, c.data.j=c.data.j, 
                                f=f, psi.inv.j = psi.inv.j, FtF=FtF)
        }
      } else {
        lmat     <- matrix(, nr=P, nc=0)
      }
      
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                                c.data=c.data, f=f, lmat=lmat)
    
      if(all(iter > burnin, iter %% thinning == 0)) {
        new.iter <- ceiling((iter - burnin)/thinning)
        psi      <- 1/psi.inv
        if(sw["mu.sw"]) mu.store[,,new.iter]    <- mu  
        if(sw["f.sw"])  f.store[,,,new.iter]    <- f # cbind(f)
        if(sw["l.sw"])  load.store[,,,new.iter] <- lmat
        if(sw["p.sw"])  psi.store[,,new.iter]   <- psi
        post.mu     <-  post.mu + mu/n.store
        post.psi    <-  post.psi + psi/n.store
        Sigma       <-  tcrossprod(lmat) + diag(psi)
        post.Sigma  <-  post.Sigma + Sigma/n.store
        log.like    <-  sum(mvdnorm(data=data, mu=mu, Sigma=Sigma, P=P, log.d=T))
        bic.mcmc    <-  max(bic.mcmc, log.like, na.rm=T)
      }  
    }
    returns   <- list(mu   = if(sw["mu.sw"]) mu.store,
                      f    = if(sw["f.sw"])  f.store, 
                      load = if(sw["l.sw"])  load.store, 
                      psi  = if(sw["p.sw"])  psi.store,
                      cov.mat    = cov.emp,
                      post.mu    = post.mu,
                      post.psi   = post.psi,
                      post.Sigma = post.Sigma,
                      bic        = 2 * bic.mcmc - pen)
    return(returns[!sapply(returns, is.null)])
  }