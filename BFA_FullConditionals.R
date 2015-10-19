# Define full conditional functions
  # Means
    sim.mu      <- function(mu.sigma, N, P, psi, data, f, load, iter, ...) {
      mu.omega  <- solve(solve(mu.sigma) + N * solve(diag(psi[,iter-1])))
      mvrnorm(mu=rep(0, P), Sigma=mu.omega) +
        mu.omega %*% solve(diag(psi[,iter-1])) %*% t(t(apply(data, 2, sum)) - t(apply(f[,,iter-1], 2, sum)) %*% t(load[,,iter-1])) }
    sim.mu      <- cmpfun(sim.mu)
    
  # Scores
    sim.omega.f <- function(Q, load, psi, iter, ...) {
      solve(diag(Q) + t(load[,,iter-1]) %*% solve(diag(psi[,iter-1])) %*% load[,,iter-1]) }
    sim.omega.f <- cmpfun(sim.omega.f)
    
    sim.scores  <- function(Q, f.omega, load, psi, data, mu, i, iter, ...) {
      mvrnorm(mu=rep(0, Q), Sigma=f.omega) + 
        (f.omega %*% t(load[,,iter-1]) %*% solve(diag(psi[,iter-1]))) %*% (data[i,] - mu[,iter]) }
    sim.scores  <- cmpfun(sim.scores)
    
  # Loadings
    sim.load    <- function(l.sigma, Q, f, psi, data, mu, j, iter, ...) {
      l.omega   <- solve(solve(l.sigma) + (1/psi[j,iter-1]) * t(f[,,iter]) %*% f[,,iter])
      t(mvrnorm(mu=rep(0, Q), Sigma=l.omega) + 
          (l.omega %*% t(f[,,iter]) * (1/psi[j,iter-1])) %*% (data[,j] - mu[j,iter])) }
    sim.load    <- cmpfun(sim.load)
    
  # Uniquenesses
    sim.psi     <- function(N, psi.alpha, psi.beta, data, mu, load, f, j, iter, ...) {
      rinvgamma(1, shape=(N+psi.alpha)/2, 
                scale=(sum(data[,j] - mu[j,iter] - load[j,,iter]%*%t(f[,,iter]))^2 + psi.beta)/2) }
    sim.psi     <- cmpfun(sim.psi)