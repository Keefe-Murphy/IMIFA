################################################
### IMIFA Full Conditionals (Shrinkage Case) ###
################################################

# Means
  sim.mu        <- function(N = NULL, P = NULL, sigma.mu = NULL, psi.inv = NULL, 
                            sum.data = NULL, sum.f = NULL, lmat = NULL, ...) {
    mu.omega    <- 1/(sigma.mu + N * psi.inv)
    U.mu        <- sqrt(mu.omega)
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- U.mu * z.mu
    mu.mu       <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f))
      as.vector(mu.mu + v.mu)
  }

# Scores
  sim.scores    <- function(N = NULL, Q = NULL, lmat = NULL,
                            psi.inv = NULL, c.data = NULL, ...) {
    load.psi    <- lmat * psi.inv
    U.f         <- chol(diag(Q) + crossprod(load.psi, lmat)) 
    z.f         <- matrix(rnorm(Q * N, 0, 1), nr=Q, nc=N)
    v.f         <- backsolve(U.f, z.f)
    mu.f        <- c.data %*% (load.psi %*% chol2inv(U.f))
      mu.f + t(v.f)
  }

# Loadings
  sim.load      <- function(D.load = NULL, Q = NULL, c.data.j = NULL, 
                            f = NULL, psi.inv.j = NULL, FtF = NULL, ...) {
    U.load      <- chol(D.load + psi.inv.j * FtF)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * chol2inv(U.load) %*% crossprod(f, c.data.j)
      mu.load + v.load
  }

# Uniquenesses
  sim.psi.inv   <- function(N = NULL, P = NULL, psi.alpha = NULL, psi.beta = NULL,
                            c.data = NULL, f = NULL, lmat = NULL, ...) { 
    rate.t      <- c.data - tcrossprod(f, lmat)
    rate.t      <- colSums(rate.t * rate.t)
    rgamma(P, shape=(N + psi.alpha)/2, 
           rate=(rate.t + psi.beta)/2) 
  }

# Local Shrinkage
  sim.phi       <- function(Q = NULL, P = NULL, phi.nu = NULL, 
                            tau = NULL, load.2 = NULL, ...) {
    rate.t      <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
    matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                  rate=rate.t), nr=P, nc=Q)
  }
    
# Global Shrinkage
  sim.delta1    <- function(Q = NULL, P = NULL, alpha.d1 = NULL, delta = NULL,
                            tau = NULL, sum.term = NULL, ...) {
    rgamma(1, shape=alpha.d1 + P * Q/2, 
              rate=1 + 0.5/delta[1] * tau %*% sum.term)
  }

  sim.deltak    <- function(Q = NULL, P = NULL, k = NULL, alpha.d2 = NULL,
                            delta = NULL, tau = NULL, sum.term = NULL) {
    rgamma(1, shape=alpha.d2 + P/2 * (Q - k + 1), 
              rate=1 + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
  }

# Priors
  # Means
    sim.mu.p    <- function(sigma.mu = NULL, P = NULL, ...) {
      U.mu      <- sqrt(1/sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu
    }
    
  # Scores
    sim.f.p     <- function(Q = NULL, N = NULL, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }
    
  # Loadings
    sim.load.p  <- function(D.load = NULL, Q = NULL, ...) {
      U.load    <- sqrt(1/D.load)
      z.load    <- rnorm(Q, 0, 1)
        U.load * z.load
    }

  # Uniquenesses
    sim.psi.p   <- function(P = NULL, psi.alpha = NULL, psi.beta = NULL, ...) {
      rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }

  # Local Shrinkage
    sim.p.p     <- function(Q = NULL, P = NULL, phi.nu = NULL, ...) {
        matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P, nc=Q)
    }

  # Global Shrinkage
    sim.d.p     <- function(Q = NULL, alpha.d1 = NULL, alpha.d2 = NULL, ...) {
        delta1  <- rgamma(n=1,     shape=alpha.d1, rate=1)
        deltak  <- rgamma(n=Q - 1, shape=alpha.d2, rate=1)
          c(delta1, deltak)
    }

  # Multivariate Normal Density
    mvdnorm     <- function(data = NULL, mu = NULL, Sigma = NULL, 
                            P = NULL, log.d = T, ...) {
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