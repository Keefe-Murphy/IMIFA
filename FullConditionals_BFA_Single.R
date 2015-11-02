####################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Single Case) ###
####################################################################################

# Preamble 
  fc.pkgs <- c("MCMCpack", "compiler")
  invisible(lapply(fc.pkgs, library, ch=T))

# Means
  sim.mu      <- function(mu.sigma, psi.inv, sum.data, sum.f, load, ...) {
    mu.omega  <- diag(1/(mu.sigma + N * psi.inv))
    mu.mu     <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
      mvrnorm(mu=rep(0, P), Sigma=mu.omega) + mu.mu 
  }; sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.omega.f <- function(Q, load, psi.inv, ...) {
    f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
    f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
      return(list(A = f.omega.a, B = f.omega.b)) 
  }; sim.omega.f <- cmpfun(sim.omega.f)
  
  sim.scores  <- function(f.omega, Q, c.data.i, ...) { 
    mu.f      <- f.omega$B %*% c.data.i
      mvrnorm(mu=rep(0, Q), Sigma=f.omega$A) + mu.f 
  }; sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load    <- function(l.sigma, Q, c.data.j, f, psi.j, FtF, ...) {
    l.omega   <- solve(l.sigma + 1/psi.j * FtF)
    mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data.j
      t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + mu.load) 
  }; sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi     <- function(c.data.j, f, load.j, ...) { 
    scale.tmp <- sum(c.data.j - tcrossprod(load.j, f))
      rinvgamma(1, shape=(N + psi.alpha)/2, 
                scale=(scale.tmp * scale.tmp + psi.beta)/2) 
  }; sim.psi     <- cmpfun(sim.psi)

# Rue & Herd Trick!!!
#
# # scores
# L.f     <- chol(f.omega$A)
# z       <- rnorm(Q, 0, 1)
# v.f     <- L.f %*% z
# mu.f    <- f.omega$B %*% (data.i - mu)
# mu.f + v.f
# 
# # loadings
# L.load  <- chol(l.omega)
# z       <- rnorm(Q, 0, 1)
# b.load  <- 1/psi.j * crossprod(f, data.j - mu.j)
# v.load  <- L.load/b.load
# m.load  <- t(L.load)/v.load
# y.load  <- t(L.load)/z
# t(y.load + m.load)
# 
# L.load  <- chol(l.omega)
# z       <- rnorm(Q, 0, 1)
# v.load  <- L.load %*% z
# mu.load <- 1/psi.j * tcrossprod(l.omega, f) %*% (data.j - mu.j)
# t(mu.load + v.load)
#
# means
# L.mu  <- chol(mu.omega)
# z     <- rnorm(P, 0, 1)
# v.mu  <- crossprod(L.mu, z)
# mu.mu <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(colSums(data)) - t(load %*% colSums(f))))
# mu.mu + v.mu