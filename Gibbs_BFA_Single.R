################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
# Preamble
  source(paste(dataDirectory, "/IMIFA-GIT/FullConditionals_BFA_Single.R", sep=""))
  if(any(range.Q) >= P) stop ("Number of factors must be less than the number of variables")

# Gibbs Sampler Function
  gibbs.single   <- function(data=data, n.iters=50000, Q=2, 
                             burnin=(n.iters/5) - 1, thin=2, scaling=T, ...) {
        
  # Remove non-numeric columns
    data       <- data[sapply(data,is.numeric)]
  
  # Centre the data (optional)
    if (scaling) {data <- scale(data, center=T, scale=F)} else  {data <- as.matrix(data)}
  
  # Define & initialise variables
    store      <- ceiling((n.iters - burnin)/thin)
    mu.store   <- matrix(NA, nr=P, nc=store);    rownames(mu.store)   <- colnames(data) 
    f.store    <- array(NA, dim=c(N, Q, store)); colnames(f.store)    <- paste("Factor",1:Q)
    load.store <- array(NA, dim=c(P, Q, store)); rownames(load.store) <- colnames(data); colnames(load.store) <- paste("Factor",1:Q)
    psi.store  <- matrix(NA, nr=P, nc=store);    rownames(psi.store)  <- colnames(data)
    mu         <- mu.store[,1]    <- mvrnorm(mu=rep(0, P), Sigma=sigma.mu * diag(P))             
    f          <- f.store[,,1]    <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))         
    load       <- load.store[,,1] <- mvrnorm(n=P, mu=rep(0, Q), Sigma=sigma.l * diag(Q))         
    psi        <- psi.store[,1]   <- rinvgamma(n=P, shape=psi.alpha/2, scale=psi.beta/2) 
    mu.sigma   <- 1/sigma.mu
    l.sigma    <- 1/sigma.l * diag(Q)
    psi.inv    <- 1/psi
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(iter < n.iters/5 && iter %% (n.iters/100) == 0) {
        cat(paste0("Iteration: ", iter, "\n"))
        } else if (iter %% (n.iters/10) == 0) {
        cat(paste0("Iteration: ", iter, "\n"))
      }
      
      sum.data   <- colSums(data)
      sum.f      <- colSums(f)
      mu         <- sim.mu(mu.sigma, psi.inv, sum.data, sum.f, load)
      f.omega    <- sim.omega.f(Q, load, psi.inv)
            
      for (i in 1:N) {
        c.data.i <- data[i,] - mu
        f[i,]    <- sim.scores(f.omega, Q, c.data.i)
      }
      FtF        <- crossprod(f)
        
      for (j in 1:P) {
        psi.j    <- psi[j]
        c.data.j <- data[,j] - mu[j]
        load[j,] <- load.j <- sim.load(l.sigma, Q, c.data.j, f, psi.j, FtF)
        psi[j]   <- sim.psi(c.data.j, f, load.j)
      } 
        psi.inv  <- 1/psi
      
      if(iter > burnin && iter %% thin == 0) {
        new.iter <- ceiling((iter-burnin)/thin)
        mu.store[,new.iter]    <- mu  
        f.store[,,new.iter]    <- f
        load.store[,,new.iter] <- load
        psi.store[,new.iter]   <- psi
      }  
    }
  return(list(mu      = mu.store,
              f       = f.store, 
              load    = load.store, 
              psi     = psi.store,
              n.store = store))
  }; gibbs.single    <- cmpfun(gibbs.single)