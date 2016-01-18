###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.data      <- function(P=25, Q=2, N=1500,
                          mean.mu=1, mean.load=0, shape.psi=1, rate.psi=0.3) {
   # N.old    <- N
   # G        <- 3
   # rand     <- sample(1:(N.old - 1), G, replace=F)
   # N.grp    <- round(rand/sum(rand) * N.old); N <- sum(N.grp)
   # i in 1:N.grp[1]
   # i in (N.grp[1] + 1):(N.grp[1] + N.grp[2])
   # i in (N.grp[1] + N.grp[2] + 1):N
  
    SimData   <- matrix(0, nr=N, nc=P)
  
  # Simulate true parameter values
    mu.true   <- mvrnorm(mu=rep(mean.mu, P), Sigma=diag(P))
    load.true <- mvrnorm(n=P, mu=rep(mean.load, Q), Sigma=diag(Q))        
    psi.true  <- 1/rgamma(n=P, shape=shape.psi, rate=rate.psi)                            
    names(mu.true)      <- c(1:P)
    colnames(load.true) <- paste("Factor", 1:Q)
    rownames(load.true) <- c(1:P)
    names(psi.true)     <- c(1:P)
  
  # Simulate data
    SimData   <- mvrnorm(n=N, mu=mu.true, Sigma=tcrossprod(load.true, load.true) + diag(psi.true))
  
  # Post-process data
    SimData   <- SimData[sample(1:N, N, replace=F),]
    rownames(SimData)   <- c(1:N)
    colnames(SimData)   <- c(1:P)
    SimData   <- as.data.frame(SimData)
    attr(SimData, "mu.true")   <- mu.true
    attr(SimData, "load.true") <- load.true
    attr(SimData, "psi.true")  <- psi.true
    return(SimData)
  }