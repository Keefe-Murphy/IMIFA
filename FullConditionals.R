###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    sim.mu      <- function(N, P, mu.sigma, psi.inv, sum.data, sum.f, lmat, mu.zero = 0) {
      mu.omega  <- 1/(mu.sigma + N * psi.inv)
      U.mu      <- sqrt(mu.omega)
      z.mu      <- rnorm(P, 0, 1)
      v.mu      <- U.mu * z.mu
      mu.mu     <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f) + mu.sigma * mu.zero)
        as.vector(mu.mu + v.mu)
    }
  
  # Scores
    sim.score   <- function(N, Q, lmat, psi.inv, c.data) {
      load.psi  <- lmat * psi.inv
      U.f       <- chol(diag(Q) + crossprod(load.psi, lmat))
      z.f       <- matrix(rnorm(Q * N, 0, 1), nr=Q, nc=N)
      v.f       <- backsolve(U.f, z.f)
      mu.f      <- c.data %*% (load.psi %*% chol2inv(U.f))
        mu.f + t(v.f)
    }
  
  # Loadings
    sim.load    <- function(l.sigma, Q, c.data, P, f, phi, tau, psi.inv, FtF, shrink = T) {
      if(shrink) {
        U.load  <- lapply(seq_len(P), function(j) chol((phi[j,] * tau * diag(Q)) + psi.inv[j] * FtF))
      } else     {
        U.load  <- lapply(seq_len(P), function(j) chol(l.sigma + psi.inv[j] * FtF))
      }
      z.load    <- lapply(seq_len(P), function(j) rnorm(Q, 0, 1))
      v.load    <- do.call(rbind, lapply(seq_len(P), function(j) backsolve(U.load[[j]], z.load[[j]])))
      mu.load   <- do.call(cbind, lapply(seq_len(P), function(j) psi.inv[j] * chol2inv(U.load[[j]]) %*% crossprod(f, c.data[,j])))
        t(mu.load) + v.load
    }

  # Uniquenesses
    sim.psi.i   <- function(N, P, psi.alpha, psi.beta, c.data, f, lmat) { 
      rate.t    <- c.data - tcrossprod(f, lmat)
      rate.t    <- colSums(rate.t * rate.t)
        rgamma(P, shape=N/2 + psi.alpha, 
               rate=rate.t/2 + psi.beta) 
    }

  # Local Shrinkage
    sim.phi     <- function(Q, P, phi.nu, tau, load.2) {
      rate.t    <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
        matrix(rgamma(P * Q, shape=1/2 + phi.nu, 
                      rate=rate.t), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta1  <- function(Q, P, alpha.d1, delta, beta.d1, tau, sum.term) {
        rgamma(1, shape=alpha.d1 + P * Q/2, 
               rate=beta.d1 + 0.5/delta[1] * tau %*% sum.term)
    }
    
    sim.deltak  <- function(Q, P, k, alpha.dk, beta.dk, delta, tau, sum.term) {
        rgamma(1, shape=alpha.dk + P/2 * (Q - k + 1), 
               rate=beta.dk + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
    }

  # Mixing Proportions
    sim.pi      <- function(pi.alpha, nn = 0) {
        rdirichlet(1, pi.alpha + nn)
    }
  
  # Cluster Labels
    sim.z       <- function(data, mu, Sigma, G, pi.prop) {
      log.numer <- do.call(cbind, lapply(seq_len(G), function(g) mvdnorm(data, mu[,g], Sigma[[g]], log.d=T) + log(pi.prop[,g])))
      log.denom <- apply(log.numer, 1, logSumExp)
      pz        <- exp(sweep(log.numer, 1, log.denom, FUN="-"))
      pz[rowSums(pz > 0) == 0,] <- rep(1/G, G)
      z         <- unlist(lapply(seq_along(log.denom), function(i) which(rmultinom(1, size=1, prob=pz[i,]) == 1)))
     #z         <- unlist(lapply(seq_along(log.denom), function(i) sample(seq_len(G), size=1, prob=pz[i,])))
        return(list(z = z, log.likes = log.denom))
    }

# Priors

  # Means
    sim.mu.p    <- function(P, sigma.mu, mu.zero = 0) {
      U.mu      <- sqrt(sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu + mu.zero
    }
  
  # Scores
    sim.f.p     <- function(Q, N) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }
  
  # Loadings
    sim.load.p  <- function(Q, P, sigma.l, phi, tau, shrink = T) {
      if(shrink) {
        U.load  <- lapply(seq_len(P), function(j) sqrt(1/(phi[j,] * tau)))
        z.load  <- lapply(seq_len(P), function(j) rnorm(Q, 0, 1))
          do.call(rbind, lapply(seq_len(P), function(j) U.load[[j]] * z.load[[j]]))
      } else     {
        U.load  <- sqrt(sigma.l)
        z.load  <- matrix(rnorm(P * Q, 0, 1), nr=P, nc=Q)
          z.load * U.load
      }
    }
  
  # Uniquenesses
    sim.psi.ip  <- function(P, psi.alpha, psi.beta) {
        rgamma(n=P, shape=psi.alpha, rate=psi.beta) 
    }

  # Local Shrinkage
    sim.phi.p   <- function(Q, P, phi.nu) {
        matrix(rgamma(n=P * Q, shape=phi.nu, rate=phi.nu), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta.p <- function(Q, alpha.d1, alpha.dk, beta.d1, beta.dk) {
      delta1    <- rgamma(n=1,     shape=alpha.d1, rate=beta.d1)
      deltak    <- rgamma(n=Q - 1, shape=alpha.dk, rate=beta.dk)
        c(delta1, deltak)
    }

  # Cluster Labels
    sim.z.p     <- function(N, prob.z) {
      ind.mat   <- rmultinom(N, size=1, prob=prob.z)
      labs      <- which(ind.mat != 0, arr.ind=T)[,1]
        factor(labs, levels=seq_along(prob.z))
    }

# Other Functions

  # Multivariate Normal Density
    mvdnorm     <- function(data, mu, Sigma, log.d = T) {
      P         <- length(mu)
      if(all(Sigma[!diag(P)] == 0)) {
        U.Sig   <- sqrt(Sigma)
      } else {
        U.Sig   <- chol(Sigma)
      }
      rooti     <- backsolve(U.Sig, diag(P))
      quad      <- crossprod(rooti, t(data) - mu)
      quads     <- colSums(quad * quad)
      log.dens  <- - P/2 * log(2 * pi) + sum(log(diag(rooti))) - 0.5 * quads
        if(isTRUE(log.d)) {
          log.dens
        } else {
          exp(log.dens)
        }
    }