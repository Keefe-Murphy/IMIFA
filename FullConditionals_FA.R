#############################################
### IMIFA Full Conditionals (Single Case) ###
#############################################
                         
# Means
  sim.mu        <- function(N, P, sigma.mu, psi.inv, sum.data, sum.f, lmat, ...) {
    mu.omega    <- 1/(sigma.mu + N * psi.inv)
    U.mu        <- sqrt(mu.omega)
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- U.mu * z.mu
    mu.mu       <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f))
      as.vector(mu.mu + v.mu)
  }

# Scores
  sim.scores    <- function(N, Q, lmat, psi.inv, c.data, ...) {
    load.psi    <- lmat * psi.inv
    f.omega.a   <- diag(Q) + crossprod(load.psi, lmat)
    U.f         <- chol(f.omega.a)
    f.omega.b   <- load.psi %*% chol2inv(U.f)
    z.f         <- matrix(rnorm(Q * N, 0, 1), nr=Q, ncol=N)
    v.f         <- backsolve(U.f, z.f)
    mu.f        <- c.data %*% f.omega.b
      mu.f + t(v.f)
  }

# Loadings
  sim.load      <- function(l.sigma, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega     <- l.sigma + psi.inv.j * FtF
    U.load      <- chol(l.omega)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * chol2inv(U.load) %*% crossprod(f, c.data.j)
      mu.load + v.load
  }

# Uniquenesses
  sim.psi.inv   <- function(N, P, psi.alpha, psi.beta, c.data, f, lmat, ...) { 
    rate.t      <- c.data - tcrossprod(f, lmat)
    rate.t      <- colSums(rate.t * rate.t)
      rgamma(P, shape=(N + psi.alpha)/2, 
                rate=(rate.t + psi.beta)/2) 
  }

# Priors
  # Means
    sim.mu.p    <- function(P, sigma.mu, ...) {
      U.mu      <- sqrt(1/sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu
    }
  
  # Scores
    sim.f.p     <- function(Q, N, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, ncol=Q)
    }

  # Loadings
    sim.l.p     <- function(Q, P, sigma.l, ...) {
      U.l       <- sqrt(1/sigma.l)
      z.l       <- matrix(rnorm(P * Q, 0, 1), nr=P, ncol=Q)
        z.l * U.l
    }

  # Uniquenesses
    sim.pi.p    <- function(P, psi.alpha, psi.beta, ...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }