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
  }; sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.omega.f <- function(Q, load, psi.inv, ...) {
    f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
    f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
    U.f       <- chol(f.omega.a)
      return(list(B = f.omega.b, U.f = U.f))
  }; sim.omega.f <- cmpfun(sim.omega.f)
  
  sim.scores  <- function(f.omega, Q, c.data.i, ...) { 
    z.f       <- rnorm(Q, 0, 1)
    v.f       <- crossprod(f.omega$U.f, z.f)
    mu.f      <- f.omega$B %*% c.data.i
    mu.f + v.f
  }; sim.scores  <- cmpfun(sim.scores)

##system.time(for(i in 1:10000){ #algorithm 2.3
#f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
#f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
#U.f       <- chol(f.omega.a)
#z.f       <- rnorm(Q, 0, 1)
#v.f       <- crossprod(f.omega$U.f, z.f)
#mu.f      <- f.omega.b %*% c.data.i
#mu.f + v.f
##})

##system.time(for(i in 1:10000){ # algorithm 2.4
#f.omega.a <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
#f.omega.b <- chol(f.omega.a)
#f.omega.c <- tcrossprod(chol2inv(f.omega.b), load) %*% diag(psi.inv)
#U.f       <- f.omega.b
##z.f       <- rnorm(Q, 0, 1)
#v.f       <- backsolve(U.f, z.f)
#system.time(for(i in 1:10000){mu.f      <- f.omega.c %*% c.data.i
#mu.f + v.f
##})

## block update
#f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
#f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
#U.f       <- chol(f.omega.a)
#z.f       <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
#v.f       <- tcrossprod(z.f, U.f)
#c.data    <- sweep(data, 2, mu) # in gibbs function
#c.data.j  <- c.data[,j]
#mu.f      <- tcrossprod(c.data, f.omega.b)
#mu.f + v.f

# Loadings
  sim.load    <- function(l.sigma, Q, c.data.j, f, psi.j, FtF, ...) {
    l.omega   <- solve(l.sigma + 1/psi.j * FtF)
    U.load    <- chol(l.omega)
    z.load    <- rnorm(Q, 0, 1)
    v.load    <- crossprod(U.load, z.load)
    mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data.j
    t(mu.load + v.load)
  }; sim.load    <- cmpfun(sim.load)

# #system.time(for(i in 1:10000){
# l.omega   <- solve(l.sigma + 1/psi.j * FtF)
# U.load    <- chol(l.omega)
# z.load    <- rnorm(Q, 0, 1)
# v.load    <- crossprod(U.load, z.load)
# mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data.j
# t(mu.load + v.load)
# #})

# #system.time(for(i in 1:10000){
# l.omega <- l.sigma + 1/psi.j * FtF
# U.load  <- chol(l.omega)
# z.load  <- rnorm(Q, 0, 1)
# v.load  <- backsolve(U.load, z.load)
# mu.load <- 1/psi.j * tcrossprod(chol2inv(U.load), f) %*% c.data.j
# t(mu.load + v.load)
# #})

# Uniquenesses
  sim.psi     <- function(c.data.j, f, load.j, ...) { 
    scale.tmp <- sum(c.data.j - tcrossprod(load.j, f))
      rinvgamma(1, shape=(N + psi.alpha)/2, 
                scale=(scale.tmp * scale.tmp + psi.beta)/2) 
  }; sim.psi     <- cmpfun(sim.psi)