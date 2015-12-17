#######################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Shrinkage Case) ###
#######################################################################################

# Preamble 
  pkgs <- c(if(exists("pkgs")) pkgs,"MCMCpack", "compiler")
  invisible(lapply(pkgs, library, ch=T))

# Means
  sim.mu       <- function(mu.sigma, psi.inv, sum.data, sum.f, load, ...) {
    mu.omega   <- diag(1/(mu.sigma + N * psi.inv))
    U.mu       <- chol(mu.omega)
    z.mu       <- rnorm(P, 0, 1)
    v.mu       <- crossprod(U.mu, z.mu)
    mu.mu      <- crossprod(crossprod(mu.omega, diag(psi.inv)), sum.data - load %*% sum.f)
      mu.mu + v.mu
}; sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.scores   <- function(Q, load, psi.inv, c.data, ...) {
    f.omega.a  <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
    U.f        <- chol(f.omega.a)
    f.omega.b  <- tcrossprod(chol2inv(U.f), load) %*% diag(psi.inv)
    z.f        <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    v.f        <- solve(U.f, t(z.f))
    mu.f       <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
}; sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load     <- function(D.load, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega    <- diag(D.load) + psi.inv.j * FtF
    U.load     <- chol(l.omega)
    z.load     <- rnorm(Q, 0, 1)
    v.load     <- backsolve(U.load, z.load)
    mu.load    <- psi.inv.j * tcrossprod(chol2inv(U.load), f) %*% c.data.j
      t(mu.load + v.load)
}; sim.load    <- cmpfun(sim.load)

# Local Shrinkage
  sim.phi      <- function(phi.nu, Q, tau, load.2, ...) {
    rate.t     <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
    matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                  rate= rate.t), nr=P, nc=Q)
}; sim.phi     <- cmpfun(sim.phi)    
    
# Global Shrinkage
  sim.delta1   <- function(delta.a1, Q, delta, tau, sum.term, ...) {
    rgamma(1, shape=delta.a1 + (P * Q)/2, 
              rate=1 + 0.5 * (1/delta[1]) * sum(tau %*% sum.term))
}; sim.delta1  <- cmpfun(sim.delta1)
  
  sim.deltak   <- function(delta.a2, Q, k, delta, tau, sum.term) {
    rgamma(1, shape=delta.a2 + P/2 * (Q - k + 1), 
              rate=1 + 0.5 * (1/delta[k]) * sum(tau[k:Q] %*% sum.term[k:Q]))
}; sim.deltak  <- cmpfun(sim.deltak)

# Uniquenesses
  sim.psi.inv  <- function(c.data, f, load, ...) { 
    rate.t     <- c.data - tcrossprod(f, load)
    rate.t     <- colSums(rate.t * rate.t)
      rgamma(P, shape=(N + psi.alpha)/2, 
                rate=(rate.t + psi.beta)/2) 
}; sim.psi.inv <- cmpfun(sim.psi.inv)