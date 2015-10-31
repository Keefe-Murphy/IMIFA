################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Single Case) ###
################################################################
  
gibbs.single   <- function(data=data, n.iters=50000, Q=2, 
                         sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5, 
                         burnin=(n.iters/5) - 1, thin=2, scaling=T, ...) {
  
  # Remove non-numeric columns
    data       <- data[sapply(data,is.numeric)]
  
  # Centre the data (optional)
    if (scaling) {data <- scale(data, center=T, scale=F)} else  {data <- as.matrix(data)}
  
  # Define & initialise variables
    N          <- nrow(data)
    P          <- ncol(data)
    if (Q>=P) stop ("Number of factors must be less than the number of variables")
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
    l.sigma    <- 1/sigma.l
    psi.inv    <- 1/psi
  
  # Iterate
    for(iter in 2:n.iters) { 
      if(iter %% (n.iters/100) == 0) cat(paste0("Iteration: ", iter, "\n"))
      
      mu         <- sim.mu(mu.sigma, N, P, psi.inv, data, f, load)
      f.omega    <- sim.omega.f(Q, load, psi.inv)
      
      for (i in 1:N) {
        data.i   <- data[i,]
        f[i,]    <- sim.scores(f.omega, Q, data.i, mu)
      }
      
      for (j in 1:P) {
        psi.j    <- psi[j]
        data.j   <- data[,j]
        mu.j     <- mu[j]
        load[j,] <- load.j <- sim.load(l.sigma, Q, f, psi.j, data.j, mu.j)
        psi[j]   <- sim.psi(N, psi.alpha, psi.beta, data.j, mu.j, load.j, f)
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
}; gibbs.single   <- cmpfun(gibbs.single)