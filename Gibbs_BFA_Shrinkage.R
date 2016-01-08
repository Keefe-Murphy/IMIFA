###################################################################
### Gibbs Sampler for Bayesian Factor Analysis (Shrinkage Case) ###
###################################################################
  
# Preamble
  source(paste(dataDirectory, "/IMIFA-GIT/FullConditionals_BFA_Shrinkage.R", sep=""))

# Gibbs Sampler Function
  gibbs.shrink <- function(data=data, n.iters=50000, Q=min(round(5 * log(P)), P),
                           burnin=n.iters/5 - 1, thin=2, 
                           centering=T, scaling=T, print=T, 
                           adapt=T, b0=0.1, b1=0.00005, prop=3/4,
                           epsilon=ifelse(centering, 0.1, 0.01), ...) {
    
  # Warning(s)
    if(Q > P)     stop("Number of factors must be less than the number of variables")
    
  # Remove non-numeric columns & (optionally) Center/Scale the data
    data       <- data[sapply(data,is.numeric)]
    data       <- scale(data, center=centering, scale=scaling)
  
  # Define & initialise variables
    n.store    <- ceiling((n.iters - burnin)/thin)
    mu.store   <- matrix(0, nr=P, nc=n.store);    rownames(mu.store)   <- colnames(data) 
    f.store    <- array(0, dim=c(N, Q, n.store)); colnames(f.store)    <- paste("Factor", 1:Q)
    load.store <- array(0, dim=c(P, Q, n.store)); rownames(load.store) <- colnames(data); colnames(load.store) <- paste("Factor", 1:Q)
    psi.store  <- matrix(0, nr=P, nc=n.store);    rownames(psi.store)  <- colnames(data)
    Q.store    <- rep(0, n.store);                Q.store[1]           <- Q 
    
    mu         <- mvrnorm(mu=rep(0, P), Sigma=sigma.mu * diag(P))             
    f          <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))         
    phi        <- matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P)
    delta      <- c(rgamma(n=1, shape=delta.a1, rate=1), rgamma(n=Q-1, shape=delta.a2, rate=1))
    tau        <- cumprod(delta)
    load       <- matrix(0, nr=P, nc=Q)
    for(j in 1:P) {
      load[j,] <- mvrnorm(n=1, mu=rep(0, Q), Sigma=diag(1/(phi[j,] * tau)))         
    }
    psi.inv    <- rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    mu.sigma   <- 1/sigma.mu
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
        sum.f        <- colSums(f)
        mu           <- sim.mu(mu.sigma, psi.inv, sum.data, sum.f, load)
        
      # Scores
        c.data       <- sweep(data, 2, mu, FUN="-")
        f            <- sim.scores(Q, load, psi.inv, c.data)
                        
      # Loadings
        FtF          <- crossprod(f)
        for (j in 1:P) {
          psi.inv.j  <- psi.inv[j]
          c.data.j   <- c.data[,j]
          D.load     <- phi[j,] * tau
          load[j,]   <- sim.load(D.load, Q, c.data.j, f, psi.inv.j, FtF)
        } 
      
      # Uniquenesses
        psi.inv      <- sim.psi.inv(c.data, f, load)
      
      # Local Shrinkage
        load.2       <- load * load
        phi          <- sim.phi(phi.nu, Q, tau, load.2)
          
      # Global Shrinkage
        sum.term     <- diag(t(phi) %*% (load * load))
        delta[1]     <- sim.delta1(delta.a1, Q, delta, tau, sum.term)
        tau          <- cumprod(delta)
        for(k in 2:Q) { 
          delta[k]   <- sim.deltak(delta.a2, Q, k, delta, tau, sum.term)
          tau        <- cumprod(delta)      
        }
      
      # Adaptation  
        if(adapt && iter > burnin) {      
          prob       <- 1/exp(b0 + b1 * pmax(iter - burnin, 0))
          unif       <- runif(n=1, min=0, max=1)
          lind       <- colSums(abs(load) < epsilon) / P
          colvec     <- lind >= prop
          numred     <- sum(colvec)
          
          if(unif    <  prob) { # check whether to adapt or not
            if(Q < P && numred == 0) { # simulate extra columns from priors
              Q      <- Q + 1
              f      <- cbind(f, rnorm(n=N, mean=0, sd=1))         
              phi    <- cbind(phi, rgamma(n=P, shape=phi.nu/2, rate=phi.nu/2))
              delta  <- c(delta, rgamma(n=1, shape=delta.a2, rate=1))
              tau    <- cumprod(delta)
              load   <- cbind(load, rnorm(n=P, mean=0, sd=sqrt(1/phi[,Q] * 1/tau[Q])))
            } else if(numred > 0) { # remove redundant columns
              nonred <- which(colvec == 0)
              Q      <- max(Q - numred, 1)
              f      <- f[,nonred]
              phi    <- phi[,nonred]
              delta  <- delta[nonred]
              tau    <- cumprod(delta)
              load   <- load[,nonred]
            }
          }
        } 
      
      if(iter > burnin && iter %% thin == 0) {
        new.iter     <- ceiling((iter - burnin)/thin)
        mu.store[,new.iter]       <- mu  
        f.store[,1:Q,new.iter]    <- f
        load.store[,1:Q,new.iter] <- load
        psi.store[,new.iter]      <- 1/psi.inv  
        Q.store[new.iter]         <- Q
      }
    }
  return(list(mu      = mu.store,
              f       = f.store, 
              load    = load.store, 
              psi     = psi.store,
              n.store = n.store,
              Q.store = Q.store))
  }; gibbs.shrink    <- cmpfun(gibbs.shrink)