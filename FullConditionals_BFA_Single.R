####################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Single Case) ###
####################################################################################

# Preamble 
  fc.pkgs <- c("MCMCpack", "compiler")
  invisible(lapply(fc.pkgs, library, ch=T))

# Means
  sim.mu      <- function(mu.sigma, psi.inv, sum.data, sum.f, load, ...) {
    mu.omega  <- diag(1/(mu.sigma + N * psi.inv))
    U.mu      <- chol(mu.omega)
    z.mu      <- rnorm(P, 0, 1)
    v.mu      <- crossprod(U.mu, z.mu)
    mu.mu     <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
      mu.mu + v.mu
}; sim.mu     <- cmpfun(sim.mu)

# Scores
  sim.scores  <- function(Q, load, psi.inv, c.data, ...) {
    f.omega.a <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
    U.f       <- chol(f.omega.a)
    f.omega.b <- tcrossprod(chol2inv(U.f), load) %*% diag(psi.inv)
    z.f       <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    v.f       <- solve(U.f, t(z.f))
    mu.f      <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
}; sim.scores <- cmpfun(sim.scores)

# Loadings
  sim.load    <- function(l.sigma, Q, c.data.j, f, psi.j, FtF, ...) {
    l.omega   <- l.sigma + 1/psi.j * FtF
    U.load    <- chol(l.omega)
    z.load    <- rnorm(Q, 0, 1)
    v.load    <- backsolve(U.load, z.load)
    mu.load   <- 1/psi.j * tcrossprod(chol2inv(U.load), f) %*% c.data.j
      t(mu.load + v.load)
}; sim.load   <- cmpfun(sim.load)

# Uniquenesses
  sim.psi     <- function(c.data.j, f, load.j, ...) { 
    scale.tmp <- sum(c.data.j - tcrossprod(load.j, f))
      rinvgamma(1, shape=(N + psi.alpha)/2, 
                scale=(scale.tmp * scale.tmp + psi.beta)/2) 
}; sim.psi    <- cmpfun(sim.psi)