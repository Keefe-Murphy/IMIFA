###############################################
### Simulate Data (Single & Shrinkage Case) ###
###############################################

# Set dimensions etc.
  P     <- 20
  Q     <- 2
  N     <- 1500
 #N.old <- 1500
 #G     <- 3
 #rand  <- sample(1:(N.old - 1), G, replace=F)
 #N.grp <- round(rand/sum(rand) * N.old); N <- sum(N.grp)
 # i in 1:N.grp[1]
 # i in (N.grp[1] + 1):(N.grp[1] + N.grp[2])
 # i in (N.grp[1] + N.grp[2] + 1):N
  
  sim.data   <- matrix(0, nr=N, nc=P)

# Simulate true parameter values
  mu.true    <- mvrnorm(mu=rep(1, P), Sigma=diag(P));             names(mu.true)      <- c(1:P)
  load.true  <- mvrnorm(n=P, mu=rep(0, Q), Sigma=diag(Q));        colnames(load.true) <- paste("Factor", 1:Q); rownames(load.true) <- c(1:P)
  psi.true   <- 1/rgamma(n=P, 1, 0.3);                            names(psi.true)     <- c(1:P)
 #f.true     <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q)); colnames(f.true)    <- paste("Factor", 1:Q); rownames(f.true)    <- c(1:N.grp[1])
 #eps.true   <- mvrnorm(n=N, mu=rep(0, P), Sigma=diag(psi.true))

# Simulate data
  sim.data   <- mvrnorm(n=N, mu=mu.true, Sigma=tcrossprod(load.true, load.true) + diag(psi.true))

# Post-process data
  sim.data   <- sim.data[sample(1:N, N, replace=F),]
  rownames(sim.data) <- c(1:N)
  colnames(sim.data) <- c(1:P)
  sim.data   <- as.data.frame(sim.data)
###