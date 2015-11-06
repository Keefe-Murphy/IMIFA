###################################
### Simulate Data (Single Case) ###
###################################

# Preamble
  sim.pkgs <- c("MCMCpack")
  invisible(lapply(sim.pkgs, library, ch=T))

# Set dimensions etc.
  P     <- 50
  Q     <- 10
  N.old <- 1500
  G     <- 3
  rand  <- sample(1:(N.old - 1), G, replace=F)
  N.grp <- round(rand/sum(rand) * N.old); N <- sum(N.grp)
  
  data  <- matrix(0, nr=N, nc=P)

# Simulate and populate each cluster
  mu.true1    <- mvrnorm(mu=rep(1, P), Sigma=diag(P));             names(mu.true1)      <- c(1:P)
  f.true1     <- mvrnorm(n=N.grp[1], mu=rep(0, Q), Sigma=diag(Q)); colnames(f.true1)    <- paste("Factor", 1:Q); rownames(f.true1)    <- c(1:N.grp[1])
  load.true1  <- mvrnorm(n=P, mu=rep(0, Q), Sigma=diag(Q));        colnames(load.true1) <- paste("Factor", 1:Q); rownames(load.true1) <- c(1:P)
  psi.true1   <- rinvgamma(n=P, 1, 1);                             names(psi.true1)     <- c(1:P)
  eps.true1   <- mvrnorm(n=N.grp[1], mu=rep(0, P), Sigma=diag(psi.true1))
  
    for (i in 1:N.grp[1]) {  
      data[i, ] <- mu.true1 + load.true1 %*% f.true1[i,] + eps.true1[i,]
    }

  mu.true2    <- mvrnorm(mu=rep(2, P), Sigma=diag(P) * 2);             names(mu.true2)      <- c(1:P)
  f.true2     <- mvrnorm(n=N.grp[2], mu=rep(2, Q), Sigma=diag(Q));     colnames(f.true2)    <- paste("Factor", 1:Q); rownames(f.true2)    <- c(1:N.grp[2])
  load.true2  <- mvrnorm(n=P, mu=rep(2, Q), Sigma=diag(Q) * 2);        colnames(load.true2) <- paste("Factor", 1:Q); rownames(load.true2) <- c(1:P)
  psi.true2   <- rinvgamma(n=P, 2, 2);                                 names(psi.true2)     <- c(1:P)
  eps.true2   <- mvrnorm(n=N.grp[2], mu=rep(0, P), Sigma=diag(psi.true2))
  
    for (i in (N.grp[1] + 1):(N.grp[1] + N.grp[2])) {  
      data[i, ] <- mu.true2 + load.true2 %*% f.true2[i-N.grp[1],] + eps.true2[i-N.grp[1],]
    }

  mu.true3    <- mvrnorm(mu=rep(3, P), Sigma=diag(P) * 3);             names(mu.true3)      <- c(1:P)
  f.true3     <- mvrnorm(n=N.grp[3], mu=rep(0, Q), Sigma=diag(Q)); colnames(f.true3)    <- paste("Factor", 1:Q); rownames(f.true3)    <- c(1:N.grp[3])
  load.true3  <- mvrnorm(n=P, mu=rep(0, Q), Sigma=diag(Q) * 3);        colnames(load.true3) <- paste("Factor", 1:Q); rownames(load.true3) <- c(1:P)
  psi.true3   <- rinvgamma(n=P, 3, 3);                             names(psi.true3)     <- c(1:P)
  eps.true3   <- mvrnorm(n=N.grp[3], mu=rep(0, P), Sigma=diag(psi.true3))
  
    for (i in (N.grp[1] + N.grp[2] + 1):N) {  
      data[i, ] <- mu.true3 + load.true3 %*% f.true3[i - N.grp[1] - N.grp[2],] + eps.true3[i - N.grp[1] - N.grp[2],]
    }

# Post-process data
  data <- data[sample(1:N, N, replace=F),]
  rownames(data) <- c(1:N)
  colnames(data) <- c(1:P)
  data           <- as.data.frame(data)