###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

sim.imifa     <- function(N = 1000, P = 25, Q = 5, G = 1) {
   
  # N.old     <- N
  # G         <- 3
  # rand      <- sample(1:(N.old - 1), G, replace=F)
  # N.grp     <- round(rand/sum(rand) * N.old); N <- sum(N.grp)
  # i in 1:N.grp[1]
  # i in (N.grp[1] + 1):(N.grp[1] + N.grp[2])
  # i in (N.grp[1] + N.grp[2] + 1):N
  if(any(N < 0, P < 0, Q < 0, G <= 0))   stop("N, P, and Q must be strictly non-negative and G must be strictly positive")
  SimData     <- matrix(0, nr=N, nc=P)
  
# Simulate true parameter values
  U.mu        <- sqrt(max(Q, 1) * diag(P))
  z.mu        <- rnorm(P, 0, 1)
  v.mu        <- U.mu %*% z.mu
  mu.mu       <- rnorm(P, 0, max(Q, 1))
  mu.true     <- as.vector(mu.mu + v.mu)
    
  if(Q > 0) {
    U.load    <- sqrt(Q * diag(Q))
    z.load    <- matrix(rnorm(P * Q, 0, 1), nr=Q, ncol=P)
    v.load    <- t(U.load %*% z.load)
    mu.load   <- rnorm(Q, 0, Q/2)
    load.true <- mu.load + v.load
  } else {
    load.true <- matrix(, nr=P, nc=0)
  }
  
  psi.true    <- 1/rgamma(n=P, shape=4, rate=1)                            
    
  names(mu.true)      <- c(1:P)
  if(Q > 0) colnames(load.true) <- paste("Factor", 1:Q)
  rownames(load.true) <- c(1:P)
  names(psi.true)     <- c(1:P)
  
# Simulate data
  omega       <- tcrossprod(load.true, load.true) + diag(psi.true)
  if(Q > 0) {
    U.om      <- chol(omega)
  } else {
    U.om      <- sqrt(omega)
  }
  z           <- matrix(rnorm(P * N, 0, 1), nr=P, ncol=N)
  v           <- U.om %*% z
  SimData     <- t(mu.true + v)
  
# Post-process data
  SimData     <- SimData[sample(1:N, N, replace=F),]
  rownames(SimData)   <- c(1:N)
  colnames(SimData)   <- c(1:P)
  SimData     <- as.data.frame(SimData)
  attr(SimData, "mu.true")      <- mu.true
  attr(SimData, "load.true")    <- load.true
  attr(SimData, "psi.true")     <- psi.true
  class(SimData)      <- c("data.frame", "IMIFA")
  return(SimData)
}