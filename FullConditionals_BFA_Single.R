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
    z         <- rnorm(P, 0, 1)
    v.mu      <- crossprod(U.mu, z)
    mu.mu     <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
    mu.mu + v.mu
  }; sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.omega.f <- function(Q, load, psi.inv, ...) {
    f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
    f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
    U.f       <- chol(f.omega.a)
    z         <- rnorm(Q, 0, 1)
    v.f       <- crossprod(U.f, z)
      return(list(B = f.omega.b, v.f = v.f))
  }; sim.omega.f <- cmpfun(sim.omega.f)
  
  sim.scores  <- function(f.omega, Q, c.data.i, ...) { 
    mu.f      <- f.omega$B %*% c.data.i
    mu.f + f.omega$v.f
  }; sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load    <- function(l.sigma, Q, c.data.j, f, psi.j, FtF, ...) {
    l.omega   <- solve(l.sigma + 1/psi.j * FtF)
    U.load    <- chol(l.omega)
    z         <- rnorm(Q, 0, 1)
    v.load    <- crossprod(U.load, z)
    mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data.j
    t(mu.load + v.load)
  }; sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi     <- function(c.data.j, f, load.j, ...) { 
    scale.tmp <- sum(c.data.j - tcrossprod(load.j, f))
      rinvgamma(1, shape=(N + psi.alpha)/2, 
                scale=(scale.tmp * scale.tmp + psi.beta)/2) 
  }; sim.psi     <- cmpfun(sim.psi)

# Rue & Herd Trick!!!
# # Precision version
# # scores
# f.omega.a <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
# f.omega.b <- chol(f.omega.a)
# f.omega.c <- tcrossprod(chol2inv(f.omega.b), load) %*% diag(psi.inv)
# L.f       <- chol(f.omega$B)
# z         <- rnorm(Q, 0, 1)
# v.f       <- backsolve(L.f, z)
# mu.f      <- f.omega$C %*% c.data.i
# mu.f + v.f
# # loadings 
# l.omega <- l.sigma + 1/psi.j * FtF
# L.load  <- chol(l.omega)
# z       <- rnorm(Q, 0, 1)
# v.load  <- backsolve(L.load, z)
# mu.load <- 1/psi.j * tcrossprod(chol2inv(L.load), f) %*% c.data.j
# t(mu.load + v.load)
# # means
# mu.omega <- mu.sigma + N * psi.inv
# L.mu  <- chol(mu.omega)
# z     <- rnorm(P, 0, 1)
# v.mu  <- backsolve(L.mu, z)
# mu.mu <- crossprod(crossprod(chol2inv(L.mu), diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
# mu.mu + v.mu