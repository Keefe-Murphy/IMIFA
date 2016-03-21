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
    if(sw["si.sw"]) {
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
    alpha.pi     <- rep(alpha.pi, G)
    pi.prop      <- sim.pi(alpha.pi=alpha.pi)
    if(z.init == "list") {
      z          <- z.list[[G]]
    } else if(z.init == "kmeans") {
      z          <- kmeans(data, G, nstart=100)$cluster
    } else {
      z          <- sim.z.p(N=N, prob.z=pi.prop)
    }
    mu           <- sim.mu.p(P=P, sigma.mu=sigma.mu, G=G)  
    f            <- sim.f.p(Q=Q, nn=tabulate(z, nbins=G))
    lmat         <- sim.load.p(Q=Q, P=P, sigma.l=sigma.l, G=G)
    psi.inv      <- sim.psi.p(P=P, psi.alpha=psi.alpha, psi.beta=psi.beta, G=G)
    l.sigma      <- sigma.l * diag(Q)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(verbose) {
        if(all(iter < burnin, iter %% (burnin/10) == 0)) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if(iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      nn         <- tabulate(z, nbins=G)
      
    # Means
      sum.data   <- apply(data, 2, tapply, z, sum)
      sum.data   <- t(replace(sum.data, which(is.na(sum.data)), 0))
      sum.f      <- do.call(cbind, lapply(f, colSums))
      mu         <- sim.mu(N=N, P=P, sigma.mu=sigma.mu, psi.inv=psi.inv,
                           sum.data=sum.data, sum.f=sum.f, lmat=lmat, G=G)
    
    # Scores & Loadings
      c.data     <- lapply(seq_len(G), function(g) sweep(data[z == g,], 2, mu[,g]))
      if(Q > 0) {
        f        <- sim.scores(nn=nn, Q=Q, lmat=lmat, psi.inv=psi.inv, c.data=c.data)
        FtF      <- lapply(f, crossprod)
        for(j in seq_len(P)) {
          psi.inv.j <- psi.inv[j,]
          c.data.j  <- lapply(c.data, function(dat) dat[,j])
          lmat[j,,] <- sim.load(l.sigma=l.sigma, Q=Q, c.data.j=c.data.j, 
                                f=f, psi.inv.j=psi.inv.j, FtF=FtF, G=G)
        }
      } else {
        f        <- array(, dim=c(N, 0, G))
        lmat     <- array(, dim=c(P, 0, G))
      }
                  
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                                c.data=c.data, f=f, lmat=lmat, G=G)
    
    # Mixing Proportions
      pi.prop    <- sim.pi(alpha.pi=alpha.pi, nn=nn)
    
      if(all(iter > burnin, iter %% thinning == 0)) {
        new.iter <- ceiling((iter - burnin)/thinning)
        psi      <- 1/psi.inv
        if(sw["mu.sw"])             mu.store[,,new.iter]    <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,,new.iter]    <- f # cbind(f)
        if(all(sw["l.sw"], Q > 0))  load.store[,,,new.iter] <- lmat
        if(sw["si.sw"])             psi.store[,,new.iter]   <- psi
        post.mu     <-  post.mu + mu/n.store
        post.psi    <-  post.psi + psi/n.store
        Sigma       <-  tcrossprod(lmat) + diag(psi)
        post.Sigma  <-  post.Sigma + Sigma/n.store
        log.like    <-  sum(mvdnorm(data=data, mu=mu, Sigma=Sigma, P=P, log.d=T))
        bic.mcmc    <-  max(bic.mcmc, log.like, na.rm=T)
      }  
    }
    returns   <- list(mu   = if(sw["mu.sw"])             mu.store,
                      f    = if(all(sw["f.sw"], Q > 0))  f.store, 
                      load = if(all(sw["l.sw"], Q > 0))  load.store, 
                      psi  = if(sw["si.sw"])             psi.store,
                      cov.mat    = cov.emp,
                      post.mu    = post.mu,
                      post.psi   = post.psi,
                      post.Sigma = post.Sigma,
                      bic        = 2 * bic.mcmc - pen)
    return(returns[!sapply(returns, is.null)])
  }