################################################
### IMIFA Full Conditionals (Shrinkage Case) ###
################################################

# Means
  sim.mu        <- function(N, P, sigma.mu, psi.inv, sum.data, sum.f, lmat, ...) {
    mu.omega    <- diag(1/(sigma.mu + N * psi.inv))
    U.mu        <- sqrt(mu.omega)
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- U.mu %*% z.mu
    mu.mu       <- mu.omega %*% (diag(psi.inv) %*% (sum.data - lmat %*% sum.f))
      mu.mu + v.mu
};  sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.scores    <- function(N, Q, lmat, psi.inv, c.data, ...) {
    f.omega.a   <- diag(Q) + crossprod(lmat, diag(psi.inv)) %*% lmat
    U.f         <- chol(f.omega.a)
    f.omega.b   <- chol2inv(U.f) %*% crossprod(lmat, diag(psi.inv))
    z.f         <- matrix(rnorm(Q * N, 0, 1), nr=Q, ncol=N)
    v.f         <- solve(U.f, z.f)
    mu.f        <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
};  sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load      <- function(D.load, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega     <- diag(D.load) + psi.inv.j * FtF
    U.load      <- chol(l.omega)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * chol2inv(U.load) %*% crossprod(f, c.data.j)
      t(mu.load + v.load)
};  sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi.inv   <- function(N, P, psi.alpha, psi.beta, c.data, f, lmat, ...) { 
    rate.t      <- c.data - tcrossprod(f, lmat)
    rate.t      <- colSums(rate.t * rate.t)
    rgamma(P, shape=(N + psi.alpha)/2, 
           rate=(rate.t + psi.beta)/2) 
  };  sim.psi.inv <- cmpfun(sim.psi.inv)

# Local Shrinkage
  sim.phi       <- function(Q, P, phi.nu, tau, load.2, ...) {
    rate.t      <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
    matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                  rate= rate.t), nr=P, nc=Q)
};  sim.phi     <- cmpfun(sim.phi)    
    
# Global Shrinkage
  sim.delta1    <- function(Q, P, delta.a1, delta, tau, sum.term, ...) {
    rgamma(1, shape=delta.a1 + (P * Q)/2, 
              rate=1 + 0.5/delta[1] * tau %*% sum.term)
};  sim.delta1  <- cmpfun(sim.delta1)

  sim.deltak    <- function(Q, P, k, delta.a2, delta, tau, sum.term) {
    rgamma(1, shape=delta.a2 + P/2 * (Q - k + 1), 
              rate=1 + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
};  sim.deltak  <- cmpfun(sim.deltak)

# Priors
  # Means
    sim.mu.p    <- function(sigma.mu, P, ...) {
      U.mu      <- sqrt(1/sigma.mu * diag(P))
      z.mu      <- rnorm(P, 0, 1)
        U.mu %*% z.mu
    }; sim.mu.p <- cmpfun(sim.mu.p)
    
  # Scores
    sim.f.p     <- function(Q, N, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, ncol=Q)
    }; sim.f.p  <- cmpfun(sim.f.p)
    
  # Loadings
    sim.l.p     <- function(D.load, Q, ...) {
      U.l       <- sqrt(diag(1/D.load))
      z.l       <- rnorm(Q, 0, 1)
        U.l %*% z.l
    }; sim.l.p  <- cmpfun(sim.l.p)

  # Uniquenesses
    sim.pi.p    <- function(P, psi.alpha, psi.beta, ...) {
      rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }; sim.pi.p <- cmpfun(sim.pi.p)

  # Local Shrinkage
    sim.p.p     <- function(Q, P, phi.nu, ...) {
        matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P)
    }; sim.p.p  <- cmpfun(sim.p.p)

  # Global Shrinkage
    sim.d.p     <- function(Q, delta.a1, delta.a2, ...) {
        delta1  <- rgamma(n=1,     shape=delta.a1, rate=1)
        deltak  <- rgamma(n=Q - 1, shape=delta.a2, rate=1)
          c(delta1, deltak)
    }; sim.d.p  <- cmpfun(sim.d.p)