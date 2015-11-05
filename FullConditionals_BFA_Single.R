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

# old method
##system.time(for(i in 1:10000){ #  ~6 seconds
#mu.omega  <- diag(1/(mu.sigma + N * psi.inv))
#mu.mu     <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
#mvrnorm(mu=rep(0, P), Sigma=mu.omega) + mu.mu 
##})

# algorithm 2.3 (fastest!)
##system.time(for(i in 1:10000){ # < 3 seconds # gives same answer as 2.4
#mu.omega  <- diag(1/(mu.sigma + N * psi.inv))
#U.mu      <- chol(mu.omega)
#z.mu      <- rnorm(P, 0, 1)
#v.mu      <- crossprod(U.mu, z.mu)
#mu.mu     <- crossprod(crossprod(mu.omega, diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
#mu.mu + v.mu
##})

# algorithm 2.4
##system.time(for(i in 1:10000){ # >4 seconds # gives same answer as 2.3
#mu.omega  <- diag(mu.sigma + N * psi.inv)
#U.mu      <- chol(mu.omega)
#z.mu      <- rnorm(P, 0, 1)
#v.mu      <- backsolve(U.mu, z.mu)
#mu.mu     <- crossprod(crossprod(chol2inv(U.mu), diag(psi.inv)), t(t(sum.data) - t(load %*% sum.f)))
#mu.mu + v.mu
#})

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

# old method
##system.time(for(i in 1:10000){ # 5 seconds
#f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
#f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
#mu.f      <- f.omega.b %*% c.data[1,]
#mvrnorm(mu=rep(0, Q), Sigma=f.omega.a) + mu.f 
##})

#algorithm 2.3
##system.time(for(i in 1:10000){ # <3 seconds # gives slightly different answer to 2.4
#f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
#f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
#U.f       <- chol(f.omega.a)
#z.f       <- rnorm(Q, 0, 1)
#v.f       <- crossprod(U.f, z.f)
#mu.f      <- f.omega.b %*% c.data[1,]
#mu.f + v.f
##})

# algorithm 2.4 (fastest!)
##system.time(for(i in 1:10000){ # > 2 seconds #gives slightly different answer to 2.3
#f.omega.a <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
#U.f       <- chol(f.omega.a)
#f.omega.b <- tcrossprod(chol2inv(U.f), load) %*% diag(psi.inv)
#z.f       <- rnorm(Q, 0, 1)
#v.f       <- backsolve(U.f, z.f)
#mu.f      <- f.omega.b %*% c.data[1,]
#mu.f + v.f
#})

## block update (2.3) # gives same answer for first row as individual method
##system.time(for(i in 1:10000){ # >8 seconds
#f.omega.a <- solve(diag(Q) + crossprod(load, diag(psi.inv)) %*% load)
#f.omega.b <- tcrossprod(f.omega.a, load) %*% diag(psi.inv)
#U.f       <- chol(f.omega.a)
#z.fb      <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
#z.fb[1,]  <- z.f
#v.f       <- z.fb %*% U.f
#mu.f      <- tcrossprod(c.data, f.omega.b)
#mu.f + v.f
#})

## block update (2.4) # (fastest!) # gives same answer for first row as individual method
##system.time(for(i in 1:10000){ # <8 seconds 
#f.omega.a <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
#U.f       <- chol(f.omega.a)
#f.omega.b <- tcrossprod(chol2inv(U.f), load) %*% diag(psi.inv)
#z.fb     <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
#z.fb[1,] <- z.f
#v.f       <- backsolve(U.f, z.f)
#mu.f      <- tcrossprod(c.data, f.omega.b)
#t(t(mu.f) + v.f)
#})

# Loadings
  sim.load    <- function(l.sigma, Q, c.data.j, f, psi.j, FtF, ...) {
    l.omega   <- solve(l.sigma + 1/psi.j * FtF)
    U.load    <- chol(l.omega)
    z.load    <- rnorm(Q, 0, 1)
    v.load    <- crossprod(U.load, z.load)
    mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data.j
    t(mu.load + v.load)
  }; sim.load    <- cmpfun(sim.load)

# old method
##system.time(for(i in 1:10000){ # 5 seconds
#l.omega   <- solve(l.sigma + 1/psi.j * FtF)
#mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data[,j]
#t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + mu.load)
##})

# algorithm 2.3
##system.time(for(i in 1:10000){ # 3 seconds #gives slightly different answer than 2.4
#l.omega   <- solve(l.sigma + 1/psi.j * FtF)
#U.load    <- chol(l.omega)
#z.load    <- rnorm(Q, 0, 1)
#v.load    <- crossprod(U.load, z.load)
#mu.load   <- 1/psi.j * tcrossprod(l.omega, f) %*% c.data[,j]
#t(mu.load + v.load)
##})

# algorithm 2.4 (fastest!)
##system.time(for(i in 1:10000){ # >2 seconds #gives slightly different answer than 2.3
#l.omega <- l.sigma + 1/psi.j * FtF
#U.load  <- chol(l.omega)
#z.load  <- rnorm(Q, 0, 1)
#v.load  <- backsolve(U.load, z.load)
#mu.load <- 1/psi.j * tcrossprod(chol2inv(U.load), f) %*% c.data[,j]
#t(mu.load + v.load)
##})

# Uniquenesses
  sim.psi     <- function(c.data.j, f, load.j, ...) { 
    scale.tmp <- sum(c.data.j - tcrossprod(load.j, f))
      rinvgamma(1, shape=(N + psi.alpha)/2, 
                scale=(scale.tmp * scale.tmp + psi.beta)/2) 
  }; sim.psi     <- cmpfun(sim.psi)