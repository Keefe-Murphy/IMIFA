###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Gibbs Sampler Function
  gibbs.IFA     <- function(Q, data, n.iters, N, P, 
                           sigma.mu, psi.alpha, psi.beta, 
                           burnin, thinning, n.store, verbose, 
                           phi.nu, delta.a1, delta.a2, 
                           adapt, b0, b1, prop, epsilon, ...) {    
    
  # Define & initialise variables
    mu.store    <- matrix(0, nr=P, nc=n.store)
    f.store     <- array(0, dim=c(N, Q, n.store))
    load.store  <- array(0, dim=c(P, Q, n.store))
    psi.store   <- matrix(0, nr=P, nc=n.store)
    Q.store     <- c(Q, rep(0, n.store - 1))
    
    dimnames(mu.store)[[1]]   <- colnames(data)
    dimnames(f.store)[[1]]    <- rownames(data)
    dimnames(f.store)[[2]]    <- paste0("Factor ", 1:Q)
    dimnames(load.store)[[1]] <- colnames(data)
    dimnames(load.store)[[2]] <- paste0("Factor ", 1:Q)
    dimnames(psi.store)[[1]]  <- colnames(data)
    dimnames(mu.store)[[2]]   <- paste0("Iteration", 1:n.store)
    dimnames(f.store)[[3]]    <- paste0("Iteration", 1:n.store)
    dimnames(load.store)[[3]] <- paste0("Iteration", 1:n.store)
    dimnames(psi.store)[[2]]  <- paste0("Iteration", 1:n.store)
    
    sigma.mu    <- 1/sigma.mu
    mu          <- sim.mu.p(sigma.mu, P)  
    f           <- sim.f.p(Q, N)
    psi.inv     <- sim.pi.p(P, psi.alpha, psi.beta)
    phi         <- sim.p.p(Q, P, phi.nu)
    delta       <- sim.d.p(Q, delta.a1, delta.a2)
    tau         <- cumprod(delta)
    lmat        <- matrix(0, nr=P, nc=Q)
    for(j in 1:P) {
      D.load    <- phi[j,] * tau
      lmat[j,]  <- sim.l.p(D.load, Q)
    }
    
    sum.data    <- colSums(data)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(verbose) {
        if(iter <= burnin && iter %% ((burnin + 1)/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if (iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
      # Means
        sum.f        <- colSums(f)
        mu           <- sim.mu(N, P, sigma.mu, psi.inv, sum.data, sum.f, lmat)
        
      # Scores
        c.data       <- sweep(data, 2, mu, FUN="-")
        f            <- sim.scores(N, Q, lmat, psi.inv, c.data)
                        
      # Loadings
        FtF          <- crossprod(f)
        for (j in 1:P) {
          psi.inv.j  <- psi.inv[j]
          c.data.j   <- c.data[,j]
          D.load     <- phi[j,] * tau
          lmat[j,]   <- sim.load(D.load, Q, c.data.j, f, psi.inv.j, FtF)
        } 
      
      # Uniquenesses
        psi.inv      <- sim.psi.inv(N, P, psi.alpha, psi.beta, c.data, f, lmat)
      
      # Local Shrinkage
        load.2       <- lmat * lmat
        phi          <- sim.phi(Q, P, phi.nu, tau, load.2)
          
      # Global Shrinkage
        sum.term     <- diag(crossprod(phi, load.2))
        delta[1]     <- sim.delta1(Q, P, delta.a1, delta, tau, sum.term)
        tau          <- cumprod(delta)
        for(k in 2:Q) { 
          delta[k]   <- sim.deltak(Q, P, k, delta.a2, delta, tau, sum.term)
          tau        <- cumprod(delta)      
        }
      
      # Adaptation  
        if(adapt && iter > burnin) {      
          prob       <- 1/exp(b0 + b1 * pmax(iter - burnin, 0))
          unif       <- runif(n=1, min=0, max=1)
          lind       <- colSums(abs(lmat) < epsilon) / P
          colvec     <- lind >= prop
          numred     <- sum(colvec)
          
          if(unif    <  prob) {        # check whether to adapt or not
            if(Q < P && numred == 0) { # simulate extra columns from priors
              Q      <- Q + 1
              f      <- cbind(f, rnorm(n=N, mean=0, sd=1))         
              phi    <- cbind(phi, rgamma(n=P, shape=phi.nu/2, rate=phi.nu/2))
              delta  <- c(delta, rgamma(n=1, shape=delta.a2, rate=1))
              tau    <- cumprod(delta)
              lmat   <- cbind(lmat, rnorm(n=P, mean=0, sd=sqrt(1/phi[,Q] * 1/tau[Q])))
            } else if(numred > 0) {    # remove redundant columns
              nonred <- which(colvec == 0)
              Q      <- max(Q - numred, 1)
              f      <- f[,nonred]
              phi    <- phi[,nonred]
              delta  <- delta[nonred]
              tau    <- cumprod(delta)
              lmat   <- lmat[,nonred]
            }
          }
        } 
      
      if(iter > burnin && iter %% thinning == 0) {
        new.iter     <- ceiling((iter - burnin)/thinning)
        mu.store[,new.iter]       <- mu  
        f.store[,1:Q,new.iter]    <- f
        load.store[,1:Q,new.iter] <- lmat
        psi.store[,new.iter]      <- 1/psi.inv  
        Q.store[new.iter]         <- Q
      }
    }
  return(list(mu      = mu.store,
              f       = as.simple_sparse_array(f.store), 
              load    = as.simple_sparse_array(load.store), 
              psi     = psi.store,
              Q.store = Q.store))
  }