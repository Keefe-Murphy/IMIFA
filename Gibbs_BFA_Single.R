################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Preamble
  source(paste(dataDirectory, "/IMIFA-GIT/FullConditionals_BFA_Single.R", sep=""))
  if(any(range.Q) >= P) stop ("Number of factors must be less than the number of variables")

# Gibbs Sampler Function
  gibbs.single <- function(data=data, n.iters=50000, Q=2,
                           burnin=n.iters/5 - 1, thin=2, 
                           centering=T, scaling=T, print=T, ...) {
        
  # Remove non-numeric columns & (optionally) Center/Scale the data
    data       <- data[sapply(data,is.numeric)]
    data       <- scale(data, center=centering, scale=scaling)
  
  # Define & initialise variables
    n.store    <- ceiling((n.iters - burnin)/thin)
    mu.store   <- matrix(NA, nr=P, nc=n.store);    rownames(mu.store)   <- colnames(data) 
    f.store    <- array(NA, dim=c(N, Q, n.store)); colnames(f.store)    <- paste("Factor", 1:Q)
    load.store <- array(NA, dim=c(P, Q, n.store)); rownames(load.store) <- colnames(data); colnames(load.store) <- paste("Factor", 1:Q)
    psi.store  <- matrix(NA, nr=P, nc=n.store);    rownames(psi.store)  <- colnames(data)
    
    mu         <- mvrnorm(mu=rep(0, P), Sigma=sigma.mu * diag(P))             
    f          <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))         
    load       <- mvrnorm(n=P, mu=rep(0, Q), Sigma=sigma.l * diag(Q))         
    psi.inv    <- rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    mu.sigma   <- 1/sigma.mu
    l.sigma    <- 1/sigma.l * diag(Q)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(print) {
        if(iter < n.iters/10 && iter %% (n.iters/100) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if (iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
      # Means
        sum.data    <- colSums(data)
        sum.f       <- colSums(f)
        mu          <- sim.mu(mu.sigma, psi.inv, sum.data, sum.f, load)
      
      # Scores
        c.data      <- sweep(data, 2, mu, FUN="-")
        f           <- sim.scores(Q, load, psi.inv, c.data)
                
      # Loadings
        FtF         <- crossprod(f)
        for (j in 1:P) {
          psi.inv.j <- psi.inv[j]
          c.data.j  <- c.data[,j]
          load[j,]  <- sim.load(l.sigma, Q, c.data.j, f, psi.inv.j, FtF)
        }
        
      # Uniquenesses
        psi.inv     <- sim.psi.inv(c.data, f, load)
      
      if(iter > burnin && iter %% thin == 0) {
        new.iter    <- ceiling((iter-burnin)/thin)
        mu.store[,new.iter]    <- mu  
        f.store[,,new.iter]    <- f
        load.store[,,new.iter] <- load
        psi.store[,new.iter]   <- 1/psi.inv
      }  
    }
  return(list(mu      = mu.store,
              f       = f.store, 
              load    = load.store, 
              psi     = psi.store,
              n.store = n.store))
  }; gibbs.single    <- cmpfun(gibbs.single)