################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Preamble
  source(paste(dataDirectory, "/IMIFA-GIT/FullConditionals_", method, ".R", sep=""))

# Gibbs Sampler Function
  gibbs.FA     <- function(data=NULL, n.iters=NULL, Q=NULL,
                           burnin=n.iters/5 - 1, thinning=2, 
                           centering=T, scaling=T, print=T, ...) {
    
  # Warning(s)
    if(!is.element(centering, c(T, F)))  stop("Arg. must be TRUE or FALSE")
    if(!is.element(scaling,   c(T, F)))  stop("Arg. must be TRUE or FALSE")
    if(!is.element(print,     c(T, F)))  stop("Arg. must be TRUE or FALSE")
        
  # Remove non-numeric columns & (optionally) Center/Scale the data
    data       <- data[sapply(data,is.numeric)]
    data       <- scale(data, center=centering, scale=scaling)    
  
  # Define & initialise variables
    n.store    <- ceiling((n.iters - burnin)/thinning)
    mu.store   <- matrix(NA, nr=P, nc=n.store)
    f.store    <- array(NA, dim=c(N, Q, n.store))
    load.store <- array(NA, dim=c(P, Q, n.store))
    psi.store  <- matrix(NA, nr=P, nc=n.store)
  
    dimnames(mu.store)[[1]]   <- colnames(data)
    dimnames(f.store)[[1]]    <- rownames(data)
    dimnames(f.store)[[2]]    <- paste0("Factor", 1:Q)
    dimnames(load.store)[[1]] <- colnames(data)
    dimnames(load.store)[[2]] <- paste0("Factor", 1:Q)
    dimnames(psi.store)[[1]]  <- colnames(data)
    dimnames(mu.store)[[2]]   <- paste0("Iteration", 1:n.store)
    dimnames(f.store)[[3]]    <- paste0("Iteration", 1:n.store)
    dimnames(load.store)[[3]] <- paste0("Iteration", 1:n.store)
    dimnames(psi.store)[[2]]  <- paste0("Iteration", 1:n.store)
      
    mu         <- sim.mu.p()  
    f          <- sim.f.p(Q)
    lmat       <- sim.l.p(Q)
    psi.inv    <- sim.pi.p()
    l.sigma    <- 1/sigma.l * diag(Q)
    sum.data   <- colSums(data)
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(print) {
        if(iter < burnin && iter %% ((burnin + 1)/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        } else if (iter %% (n.iters/10) == 0) {
          cat(paste0("Iteration: ", iter, "\n"))
        }
      }
      
      # Means
        sum.f       <- colSums(f)
        mu          <- sim.mu(psi.inv, sum.data, sum.f, lmat)
      
      # Scores
        c.data      <- sweep(data, 2, mu, FUN="-")
        f           <- sim.scores(Q, lmat, psi.inv, c.data)
                
      # Loadings
        FtF         <- crossprod(f)
        for (j in 1:P) {
          psi.inv.j <- psi.inv[j]
          c.data.j  <- c.data[,j]
          lmat[j,]  <- sim.load(l.sigma, Q, c.data.j, f, psi.inv.j, FtF)
        }
        
      # Uniquenesses
        psi.inv     <- sim.psi.inv(c.data, f, lmat)
      
      if(iter > burnin && iter %% thinning == 0) {
        new.iter    <- ceiling((iter - burnin)/thinning)
        mu.store[,new.iter]    <- mu  
        f.store[,,new.iter]    <- f
        load.store[,,new.iter] <- lmat
        psi.store[,new.iter]   <- 1/psi.inv
      }  
    }
  return(list(mu      = mu.store,
              f       = f.store, 
              load    = load.store, 
              psi     = psi.store,
              n.store = n.store))
  }; gibbs.FA    <- cmpfun(gibbs.FA)