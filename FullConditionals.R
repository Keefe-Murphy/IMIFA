###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    sim.mu      <- function(N = NULL, P = NULL, mu.sigma = NULL, psi.inv = NULL, 
                            sum.data = NULL, sum.f = NULL, lmat = NULL, mu.zero = rep(0, P), ...) {
      mu.omega  <- 1/(mu.sigma + N * psi.inv)
      U.mu      <- sqrt(mu.omega)
      z.mu      <- rnorm(P, 0, 1)
      v.mu      <- U.mu * z.mu
      mu.mu     <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f) + mu.sigma * mu.zero)
        as.vector(mu.mu + v.mu)
    }
  
  # Scores
    sim.score   <- function(N = NULL, Q = NULL, lmat = NULL,
                            psi.inv = NULL, c.data = NULL, ...) {
      load.psi  <- lmat * psi.inv
      U.f       <- chol(diag(Q) + crossprod(load.psi, lmat))
      z.f       <- matrix(rnorm(Q * N, 0, 1), nr=Q, nc=N)
      v.f       <- backsolve(U.f, z.f)
      mu.f      <- c.data %*% (load.psi %*% chol2inv(U.f))
        mu.f + t(v.f)
    }
  
  # Loadings
    sim.load    <- function(l.sigma = NULL, Q = NULL, c.data = NULL, P = NULL, f = NULL,
                            phi = NULL, tau = NULL, psi.inv = NULL, FtF = NULL, shrink = T, ...) {
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
    sim.psi.i   <- function(N = NULL, P = NULL, psi.alpha = NULL, psi.beta = NULL, 
                            c.data = NULL, f = NULL, lmat = NULL, ...) { 
      rate.t    <- c.data - tcrossprod(f, lmat)
      rate.t    <- colSums(rate.t * rate.t)
        rgamma(P, shape=(N + psi.alpha)/2, 
               rate=(rate.t + psi.beta)/2) 
    }

  # Local Shrinkage
    sim.phi     <- function(Q = NULL, P = NULL, phi.nu = NULL, 
                            tau = NULL, load.2 = NULL, ...) {
      rate.t    <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
        matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                      rate=rate.t), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta1  <- function(Q = NULL, P = NULL, alpha.d1 = NULL, delta = NULL,
                            tau = NULL, sum.term = NULL, ...) {
        rgamma(1, shape=alpha.d1 + P * Q/2, 
               rate=1 + 0.5/delta[1] * tau %*% sum.term)
    }
    
    sim.deltak  <- function(Q = NULL, P = NULL, k = NULL, alpha.d2 = NULL,
                            delta = NULL, tau = NULL, sum.term = NULL) {
        rgamma(1, shape=alpha.d2 + P/2 * (Q - k + 1), 
               rate=1 + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
    }

  # Mixing Proportions
    sim.pi      <- function(pi.alpha = NULL, nn = 0, ...) {
        rdirichlet(1, pi.alpha + nn)
    }
  
  # Cluster Labels
    sim.z       <- function(data = NULL, mu = NULL, Sigma = NULL, 
                            G = NULL, pi.prop = NULL, ...) {
      numer     <- do.call(cbind, lapply(seq_len(G), function(g) exp(mvdnorm(data, mu[,g], Sigma[[g]], log.d=T) + log(pi.prop[,g]))))
      denomin   <- rowSums(numer)
      pz        <- sweep(numer, 1, denomin, FUN="/")
      pz[is.na(pz)]             <- 1/G
      pz[rowSums(pz > 0) == 0,] <- rep(1/G, G)
      pz[pz <= 0]               <- .Machine$double.eps
      z         <- unlist(lapply(seq_along(denomin), function(i) which(rmultinom(1, size=1, prob=pz[i,]) == 1)))
        return(list(z = z, log.likes = log(denomin)))
    }

# Priors

  # Means
    sim.mu.p    <- function(P = NULL, sigma.mu = NULL, mu.zero = rep(0, P), ...) {
      U.mu      <- sqrt(sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu + mu.zero
    }
  
  # Scores
    sim.f.p     <- function(Q = NULL, N = NULL, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }
  
  # Loadings
    sim.load.p  <- function(Q = NULL, P = NULL, sigma.l = NULL, 
                            phi = NULL, tau = NULL, shrink = T, ...) {
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
    sim.psi.ip  <- function(P = NULL, psi.alpha = NULL, psi.beta = NULL, ...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }

  # Local Shrinkage
    sim.phi.p   <- function(Q = NULL, P = NULL, phi.nu = NULL, ...) {
        matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta.p <- function(Q = NULL, alpha.d1 = NULL, alpha.d2 = NULL, ...) {
      delta1    <- rgamma(n=1,     shape=alpha.d1, rate=1)
      deltak    <- rgamma(n=Q - 1, shape=alpha.d2, rate=1)
        c(delta1, deltak)
    }

  # Cluster Labels
    sim.z.p     <- function(N = NULL, prob.z = NULL, ...) {
      ind.mat   <- rmultinom(N, size=1, prob=prob.z)
      labs      <- which(ind.mat != 0, arr.ind=T)[,1]
        factor(labs, levels=seq_along(prob.z))
    }

# Other Functions

  # Multivariate Normal Density
    mvdnorm     <- function(data = NULL, mu = NULL, Sigma = NULL, log.d = T, ...) {
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