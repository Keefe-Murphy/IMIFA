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
      mu.mu + v.mu
  }

# Scores
  sim.scores    <- function(N, Q, lmat, psi.inv, c.data, ...) {
    f.omega.a   <- diag(Q) + crossprod(lmat, diag(psi.inv)) %*% lmat
    U.f         <- chol(f.omega.a)
    f.omega.b   <- chol2inv(U.f) %*% crossprod(lmat, diag(psi.inv))
    z.f         <- matrix(rnorm(Q * N, 0, 1), nr=Q, ncol=N)
    v.f         <- solve(U.f, z.f)
    mu.f        <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
  }

# Loadings
  sim.load      <- function(l.sigma, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega     <- l.sigma + psi.inv.j * FtF
    U.load      <- chol(l.omega)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * chol2inv(U.load) %*% crossprod(f, c.data.j)
      t(mu.load + v.load)
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
      U.mu      <- sqrt(1/sigma.mu) * diag(P)
      z.mu      <- rnorm(P, 0, 1)
        U.mu %*% z.mu
    }
  
  # Scores
    sim.f.p     <- function(Q, N, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, ncol=Q)
    }

  # Loadings
    sim.l.p     <- function(Q, P, sigma.l, ...) {
      U.l       <- sqrt(1/sigma.l) * diag(Q)
      z.l       <- matrix(rnorm(P * Q, 0, 1), nr=Q, ncol=P)
        t(U.l %*% z.l)
    }

  # Uniquenesses
    sim.pi.p    <- function(P, psi.alpha, psi.beta, ...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }