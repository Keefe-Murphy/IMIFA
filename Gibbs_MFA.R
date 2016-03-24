################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Group Case) ####
################################################################
  
# Gibbs Sampler Function
  gibbs.MFA      <- function(Q = NULL, data = NULL, n.iters = NULL,
                             N = NULL, P = NULL, sigma.mu = NULL,
                             psi.alpha = NULL, psi.beta = NULL, G = NULL,
                             burnin = NULL, thinning = NULL, sw = NULL,
                             n.store = NULL, verbose = NULL, sigma.l = NULL,
                             alpha.pi = NULL, z.init = NULL, z.list = NULL,  ...) {
        
  # Define & initialise variables
    obsnames     <- rownames(data)
    varnames     <- colnames(data)
    facnames     <- paste0("Factor ", seq_len(Q))
    gnames       <- paste0("Group ", seq_len(G))
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
    post.mu      <- matrix(0, nr=P, nc=G)
    post.psi     <- matrix(0, nr=P, nc=G)
   #post.Sigma   <- array(0, dim=c(P, P, G))
    ll.max       <- - Inf
    K            <- G - 1 + G * (P * Q - 0.5 * Q * (Q - 1)) + 2 * G * P
   #cov.emp      <- cov(data)
    dimnames(post.mu)      <- dimnames(post.psi)   <- list(varnames, gnames)
   #dimnames(post.Sigma)   <- list(varnames, varnames, gnames)
   #dimnames(cov.emp)      <- dimnames(post.Sigma)
    
    sigma.mu     <- 1/sigma.mu
    sigma.l      <- 1/sigma.l
    pi.alpha     <- rep(alpha.pi, G)
    pi.prop      <- sim.pi(pi.alpha=pi.alpha)
    if(z.init == "list") {
      z          <- z.list[[G]]
    } else if(z.init == "kmeans") {
      z          <- factor(kmeans(data, G, nstart=100)$cluster, levels=seq_len(G))
    } else {
      z          <- sim.z.p(N=N, prob.z=pi.prop)
    }
    mu           <- sim.mu.p(P=P, sigma.mu=sigma.mu, G=G) 
    f            <- sim.f.p(Q=Q, N=N)
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
      if(G > 1) {
        sum.data <- t(replace(sum.data, which(is.na(sum.data)), 0))
      }
      sum.f      <- do.call(cbind, lapply(seq_len(G), function(g) colSums(f[z == g,, drop=F])))
      mu         <- sim.mu(nn=nn, P=P, sigma.mu=sigma.mu, psi.inv=psi.inv,
                           sum.data=sum.data, sum.f=sum.f, lmat=lmat, G=G)
    
    # Scores & Loadings
      c.data     <- lapply(seq_len(G), function(g) sweep(data[z == g,], 2, mu[,g], FUN="-"))
      if(Q > 0) {
        f        <- sim.scores(nn=nn, Q=Q, lmat=lmat, psi.inv=psi.inv, 
                               c.data=c.data)[as.character(seq_len(N)),, drop=F]
        FtF      <- lapply(seq_len(G), function(g) crossprod(f[z == g,, drop=F]))
        for(j in seq_len(P)) {
          psi.inv.j <- psi.inv[j,]
          c.data.j  <- lapply(c.data, function(dat) dat[,j])
          lmat[j,,] <- sim.load(l.sigma=l.sigma, Q=Q, c.data.j=c.data.j, 
                                f=f, psi.inv.j=psi.inv.j, FtF=FtF, G=G, z=z)
        }
      } else {
        f        <- matrix(, nr=N, nc=0)
        lmat     <- array(, dim=c(P, 0, G))
      }
                  
    # Uniquenesses
      psi.inv    <- sim.psi.inv(N=N, P=P, psi.alpha=psi.alpha, psi.beta=psi.beta,
                                c.data=c.data, f=f, lmat=lmat, G=G, z=z)
    
    # Mixing Proportions
      pi.prop    <- sim.pi(pi.alpha=pi.alpha, nn=nn)
    
    # Cluster Labels
      psi        <- 1/psi.inv
      Sigma      <- lapply(seq_len(G), function(g) tcrossprod(adrop(lmat[,,g, drop=F], drop=3)) + diag(psi[,g]))
      z.res      <- sim.z(data=data, mu=mu, Sigma=Sigma, 
                          N=N, G=G, P=P, pi.prop=pi.prop)
      z          <- z.res$z
    
      if(all(iter > burnin, iter %% thinning == 0)) {
        new.iter <- ceiling((iter - burnin)/thinning)
        if(sw["mu.sw"])             mu.store[,,new.iter]    <- mu  
        if(all(sw["f.sw"], Q > 0))  f.store[,,new.iter]     <- f
        if(all(sw["l.sw"], Q > 0))  load.store[,,,new.iter] <- lmat
        if(sw["psi.sw"])            psi.store[,,new.iter]   <- psi
        if(sw["pi.sw"])             pi.store[,new.iter]     <- pi.prop
                                    z.store[,new.iter]      <- z 
        post.mu     <- post.mu + mu/n.store
        post.psi    <- post.psi + psi/n.store
       #post.Sigma  <- post.Sigma + Sigma/n.store
        log.like    <- sum(z.res$log.liks)
        ll.max      <- max(ll.max, log.like, na.rm=T)
      }  
    }
    post.z    <- setNames(apply(z.store, 1, function(x) factor(which.max(tabulate(x)), levels=seq_len(G))), obsnames)
    post.pi   <- setNames(prop.table(tabulate(post.z, nbins=G)), gnames)
    returns   <- list(mu      = if(sw["mu.sw"])             mu.store,
                      f       = if(all(sw["f.sw"], Q > 0))  f.store, 
                      load    = if(all(sw["l.sw"], Q > 0))  load.store, 
                      psi     = if(sw["psi.sw"])            psi.store,
                      pi.prop = if(sw["pi.sw"])             pi.store,
                      z          = z.store,
                      post.mu    = post.mu,
                      post.psi   = post.psi,
                      post.z     = post.z,
                      post.pi    = post.pi,
                     #cov.mat    = cov.emp,
                     #post.Sigma = post.Sigma,
                      aic        = 2 * ll.max - K * 2,
                      bic        = 2 * ll.max - K * log(N))
    return(returns[!sapply(returns, is.null)])
  }