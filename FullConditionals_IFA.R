#######################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Shrinkage Case) ###
#######################################################################################

# Set hyperparameter values
  if(!exists("sigma.mu"))  assign("sigma.mu",  0.5)
  if(!exists("sigma.l"))   assign("sigma.l",   0.5)
  if(!exists("psi.alpha")) assign("psi.alpha", 2)
  if(!exists("psi.beta"))  assign("psi.beta",  0.6)
  if(!exists("phi.nu"))    assign("phi.nu",    3)
  if(!exists("delta.a1"))  assign("delta.a1",  2.1)
  if(!exists("delta.a2"))  assign("delta.a2",  12.1)
                           assign("mu.sigma",  1/sigma.mu)

# Means
  sim.mu        <- function(psi.inv, sum.data, sum.f, lmat, ...) {
    mu.omega    <- diag(1/(mu.sigma + N * psi.inv))
    U.mu        <- chol(mu.omega)
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- crossprod(U.mu, z.mu)
    mu.mu       <- crossprod(crossprod(mu.omega, diag(psi.inv)), sum.data - lmat %*% sum.f)
      mu.mu + v.mu
};  sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.scores    <- function(Q, lmat, psi.inv, c.data, ...) {
    f.omega.a   <- diag(Q) + crossprod(lmat, diag(psi.inv)) %*% lmat
    U.f         <- chol(f.omega.a)
    f.omega.b   <- tcrossprod(chol2inv(U.f), lmat) %*% diag(psi.inv)
    z.f         <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    v.f         <- solve(U.f, t(z.f))
    mu.f        <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
};  sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load      <- function(D.load, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega     <- diag(D.load) + psi.inv.j * FtF
    U.load      <- chol(l.omega)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * tcrossprod(chol2inv(U.load), f) %*% c.data.j
      t(mu.load + v.load)
};  sim.load    <- cmpfun(sim.load)

# Local Shrinkage
  sim.phi       <- function(Q, tau, load.2, ...) {
    rate.t      <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
    matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                  rate= rate.t), nr=P, nc=Q)
};  sim.phi     <- cmpfun(sim.phi)    
    
# Global Shrinkage
  sim.delta1    <- function(Q, delta, tau, sum.term, ...) {
    rgamma(1, shape=delta.a1 + (P * Q)/2, 
              rate=1 + 0.5 * (1/delta[1]) * sum(tau %*% sum.term))
};  sim.delta1  <- cmpfun(sim.delta1)
  
  sim.deltak    <- function(Q, k, delta, tau, sum.term) {
    rgamma(1, shape=delta.a2 + P/2 * (Q - k + 1), 
              rate=1 + 0.5 * (1/delta[k]) * sum(tau[k:Q] %*% sum.term[k:Q]))
};  sim.deltak  <- cmpfun(sim.deltak)

# Uniquenesses
  sim.psi.inv   <- function(c.data, f, lmat, ...) { 
    rate.t      <- c.data - tcrossprod(f, lmat)
    rate.t      <- colSums(rate.t * rate.t)
      rgamma(P, shape=(N + psi.alpha)/2, 
                rate=(rate.t + psi.beta)/2) 
};  sim.psi.inv <- cmpfun(sim.psi.inv)

# Priors
  # Means
    sim.mu.p    <- function(...) {
      mu.omega  <- sigma.mu * diag(P)
      U.mu      <- sqrt(mu.omega)
      z.mu      <- rnorm(P, 0, 1)
        as.vector(crossprod(U.mu, z.mu))
    }; sim.mu.p <- cmpfun(sim.mu.p)
    
  # Scores
    sim.f.p     <- function(Q, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, ncol=Q)
    }; sim.f.p  <- cmpfun(sim.f.p)
    
  # Loadings
    sim.l.p     <- function(D.load, Q, ...) {
      l.omega   <- diag(1/D.load)
      U.l       <- sqrt(l.omega)
      z.l       <- rnorm(Q, 0, 1)
        crossprod(U.l, z.l)
    }; sim.l.p  <- cmpfun(sim.l.p)

  # Local Shrinkage
    sim.p.p     <- function(Q, ...) {
        matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P)
    }; sim.p.p  <- cmpfun(sim.p.p)

  # Global Shrinkage
    sim.d.p     <- function(Q, ...) {
        delta1  <- rgamma(n=1, shape=delta.a1, rate=1)
        deltak  <- rgamma(n=Q - 1, shape=delta.a2, rate=1)
          c(delta1, deltak)
    }; sim.d.p  <- cmpfun(sim.d.p)

  # Uniquenesses
    sim.pi.p    <- function(...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }; sim.pi.p <- cmpfun(sim.pi.p)