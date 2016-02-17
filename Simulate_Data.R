###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.imifa     <- function(N = 1000, P = 25, Q = 5) {
   
  # N.old     <- N
  # G         <- 3
  # rand      <- sample(1:(N.old - 1), G, replace=F)
  # N.grp     <- round(rand/sum(rand) * N.old); N <- sum(N.grp)
  # i in 1:N.grp[1]
  # i in (N.grp[1] + 1):(N.grp[1] + N.grp[2])
  # i in (N.grp[1] + N.grp[2] + 1):N
  
    SimData   <- matrix(0, nr=N, nc=P)
  
  # Simulate true parameter values
    U.mu      <- sqrt(Q * diag(P))
    z.mu      <- rnorm(P, 0, 1)
    v.mu      <- U.mu %*% z.mu
    mu.mu     <- rnorm(P, 0, Q)
    mu.true   <- as.vector(mu.mu + v.mu)
    
    U.load    <- sqrt(Q * diag(Q))
    z.load    <- matrix(rnorm(P * Q, 0, 1), nr=Q, ncol=P)
    v.load    <- t(U.load %*% z.load)
    mu.load   <- rnorm(Q, 0, Q/2)
    load.true <- mu.load + v.load
    
    psi.true  <- 1/rgamma(n=P, shape=4, rate=1)                            
    
    names(mu.true)      <- c(1:P)
    colnames(load.true) <- paste("Factor", 1:Q)
    rownames(load.true) <- c(1:P)
    names(psi.true)     <- c(1:P)
  
  # Simulate data
    omega     <- tcrossprod(load.true, load.true) + diag(psi.true)
    U         <- chol(omega)
    z         <- matrix(rnorm(P * N, 0, 1), nr=P, ncol=N)
    v         <- U %*% z
    SimData   <- t(mu.true + v)
  
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