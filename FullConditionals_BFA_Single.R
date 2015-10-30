####################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Single Case) ###
####################################################################################
  
# Means
  sim.mu      <- function(mu.sigma, N, P, psi, data, f, load, ...) {
    mu.omega  <- solve(solve(mu.sigma) + N * solve(diag(psi)))
    mvrnorm(mu=rep(0, P), Sigma=mu.omega) +
      crossprod(mu.omega, solve(diag(psi))) %*% t(t(colSums(data)) - t(colSums(f)) %*% t(load)) }
  sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.omega.f <- function(Q, load, psi, ...) {
    f.omega.a <- solve(diag(Q) + crossprod(load, solve(diag(psi))) %*% load)
    f.omega.b <- tcrossprod(f.omega.a, load) %*% solve(diag(psi))
    return(list(A = f.omega.a, B = f.omega.b)) }
  sim.omega.f <- cmpfun(sim.omega.f)

  sim.scores  <- function(f.omega, Q, data.i, mu, ...) { 
    mvrnorm(mu=rep(0, Q), Sigma=f.omega$A) + 
      f.omega$B %*% (data.i - mu) }
  sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load    <- function(l.sigma, Q, f, psi.j, data.j, mu.j, ...) {
    l.omega   <- solve(solve(l.sigma) + (1/psi.j) * crossprod(f))
    t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + 
        (1/psi.j * tcrossprod(l.omega, f)) %*% (data.j - mu.j)) }
  sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi     <- function(N, psi.alpha, psi.beta, data.j, mu.j, load.j, f, ...) { 
    scale.tmp <- sum(data.j - mu.j - tcrossprod(load.j, f))
    rinvgamma(1, shape=(N + psi.alpha)/2, 
              scale=(scale.tmp * scale.tmp + psi.beta)/2) }
  sim.psi     <- cmpfun(sim.psi)