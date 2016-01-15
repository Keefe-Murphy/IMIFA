####################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Single Case) ###
####################################################################################

# Set hyperparameter values
  if(!exists("n.iters"))   assign("n.iters",   50000)
  if(!exists("N"))         assign("N",         nrow(data))
  if(!exists("P"))         assign("P",         sum(sapply(data, is.numeric)))
  if(!exists("sigma.mu"))  assign("sigma.mu",  0.5)
  if(!exists("sigma.l"))   assign("sigma.l",   0.5)
  if(!exists("psi.alpha")) assign("psi.alpha", 2)
  if(!exists("psi.beta"))  assign("psi.beta",  0.6)
                           assign("mu.sigma",  1/sigma.mu)
                           
# Means
  sim.mu        <- function(psi.inv, sum.data, sum.f, load, ...) {
    mu.omega    <- diag(1/(mu.sigma + N * psi.inv))
    U.mu        <- chol(mu.omega)
    z.mu        <- rnorm(P, 0, 1)
    v.mu        <- crossprod(U.mu, z.mu)
    mu.mu       <- crossprod(crossprod(mu.omega, diag(psi.inv)), sum.data - load %*% sum.f)
      mu.mu + v.mu
};  sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.scores    <- function(Q, load, psi.inv, c.data, ...) {
    f.omega.a   <- diag(Q) + crossprod(load, diag(psi.inv)) %*% load
    U.f         <- chol(f.omega.a)
    f.omega.b   <- tcrossprod(chol2inv(U.f), load) %*% diag(psi.inv)
    z.f         <- mvrnorm(n=N, mu=rep(0, Q), Sigma=diag(Q))
    v.f         <- solve(U.f, t(z.f))
    mu.f        <- tcrossprod(f.omega.b, c.data)
      t(mu.f + v.f)
};  sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load      <- function(l.sigma, Q, c.data.j, f, psi.inv.j, FtF, ...) {
    l.omega     <- l.sigma + psi.inv.j * FtF
    U.load      <- chol(l.omega)
    z.load      <- rnorm(Q, 0, 1)
    v.load      <- backsolve(U.load, z.load)
    mu.load     <- psi.inv.j * tcrossprod(chol2inv(U.load), f) %*% c.data.j
      t(mu.load + v.load)
};  sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi.inv   <- function(c.data, f, load, ...) { 
    rate.t      <- c.data - tcrossprod(f, load)
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
    sim.l.p     <- function(Q, ...) {
      l.omega   <- sigma.l * diag(Q)
      U.l       <- sqrt(l.omega)
      z.l       <- matrix(rnorm(P * Q, 0, 1), nr=Q, ncol=P)
        t(crossprod(U.l, z.l))
    }; sim.l.p  <- cmpfun(sim.l.p)

  # Uniquenesses
    sim.pi.p    <- function(...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }; sim.pi.p <- cmpfun(sim.pi.p)