# Bayesian Factor Analysis (Single Case)
  gibbs.single  <- function(data, n.iters=100000, Q=2, sigma.mu=0.5, sigma.l=0.5, psi.alpha=5, psi.beta=5, scaling=T, ...) {
    
  # Remove non-numeric columns
    data      <- data[sapply(data,is.numeric)]
    
  # Centre the data (optional)
    if (scaling) {data <- scale(data, center=T, scale=F)} else  {data <- as.matrix(data)}
    
  # Define & initialise variables
    N         <- nrow(data)
    P         <- ncol(data)
    if (Q>=P) stop ("Number of factors must be less than the number of variables")
    mu        <- matrix(NA, nr=P, nc=n.iters);    rownames(mu)   <- colnames(data)
    f         <- array(NA, dim=c(N, Q, n.iters)); colnames(f)    <- paste("Factor",1:Q)
    load      <- array(NA, dim=c(P, Q, n.iters)); rownames(load) <- colnames(data); colnames(load) <- paste("Factor",1:Q)
    psi       <- matrix(NA, nr=P, nc=n.iters);    rownames(psi)  <- colnames(data)
    mu.sigma  <- sigma.mu * diag(P)
    l.sigma   <- sigma.l  * diag(Q)
    mu.omega  <- matrix(NA, P, P)
    f.omega   <- matrix(NA, Q, Q)
    l.omega   <- matrix(NA, Q, Q)
    mu[,1]    <- mvrnorm(mu=rep(0, P), Sigma=mu.sigma)
    f[,,1]    <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    load[,,1] <- mvrnorm(n=P, mu=rep(0, Q), Sigma=l.sigma)
    psi[,1]   <- rinvgamma(P, shape=psi.alpha/2, scale=psi.beta/2)
    
  # Iterate
    for(iter in 2:n.iters) { 
      mu[,iter]       <- sim.mu(mu.sigma, N, P, psi, data, f, load, iter)
      f.omega         <- sim.omega.f(Q, load, psi, iter)
      for (i in 1:N) {
        f[i,,iter]    <- sim.scores(Q, f.omega, load, psi, data, mu, i, iter)
      }
      for (j in 1:P) {
        load[j,,iter] <- sim.load(l.sigma, Q, f, psi, data, mu, j, iter)
        psi[j,iter]   <- sim.psi(N, psi.alpha, psi.beta, data, mu, load, f, j, iter)
      }
    }
    return(list(mu   = mu,
                f    = f, 
                load = load, 
                psi  = psi))
  }; gibbs.single   <- cmpfun(gibbs.single)