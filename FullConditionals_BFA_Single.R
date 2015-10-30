####################################################################################
### Define full conditional functions for Bayesian Factor Analysis (Single Case) ###
####################################################################################
  
# Means
  sim.mu      <- function(mu.sigma, N, P, psi, data, f, load, ...) {
    mu.omega  <- solve(solve(mu.sigma) + N * solve(diag(psi)))
    mvrnorm(mu=rep(0, P), Sigma=mu.omega) +
      mu.omega %*% solve(diag(psi)) %*% t(t(apply(data, 2, sum)) - t(apply(f, 2, sum)) %*% t(load)) }
  sim.mu      <- cmpfun(sim.mu)

# Scores
  sim.omega.f <- function(Q, load, psi, ...) {
    f.omega.a <- solve(diag(Q) + t(load) %*% solve(diag(psi)) %*% load) 
    f.omega.b <- f.omega.a %*% t(load) %*% solve(diag(psi)) 
    return(list(A = f.omega.a, B = f.omega.b)) }
  sim.omega.f <- cmpfun(sim.omega.f)

  sim.scores  <- function(f.omega, Q, data.i, mu, ...) { 
    mvrnorm(mu=rep(0, Q), Sigma=f.omega$A) + 
      f.omega$B %*% (data.i - mu) }
  sim.scores  <- cmpfun(sim.scores)

# Loadings
  sim.load    <- function(l.sigma, Q, f, psi.j, data.j, mu.j, ...) {
    l.omega   <- solve(solve(l.sigma) + (1/psi.j) * t(f) %*% f)
    t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + 
        (l.omega %*% t(f) * (1/psi.j)) %*% (data.j - mu.j)) }
  sim.load    <- cmpfun(sim.load)

# Uniquenesses
  sim.psi     <- function(N, psi.alpha, psi.beta, data.j, mu.j, load.j, f, ...) { 
    scale.tmp <- sum(data.j - mu.j - load.j %*% t(f))
    rinvgamma(1, shape=(N + psi.alpha)/2, 
              scale=(scale.tmp * scale.tmp + psi.beta)/2) }
  sim.psi     <- cmpfun(sim.psi)