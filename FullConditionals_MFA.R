#############################################
### IMIFA Full Conditionals (Single Case) ###
#############################################
                         
# Means
  sim.mu        <- function(N = NULL, P = NULL, sigma.mu = NULL, psi.inv = NULL, 
                            sum.data = NULL, sum.f = NULL, lmat = NULL, G = NULL, ...) {
    mu.omega    <- 1/(sigma.mu + N * psi.inv)
    U.mu        <- apply(mu.omega, 2, sqrt)
    z.mu        <- matrix(rnorm(P * G, 0, 1), nr=P, nc=G)
    v.mu        <- U.mu * z.mu
    lf.prod     <- do.call(cbind, lapply(seq_len(G), function(g) lmat[,,g] %*% sum.f[,g]))
    mu.mu       <- mu.omega * (psi.inv * (sum.data - lf.prod))
      mu.mu + v.mu
  }

# Scores
  sim.scores    <- function(nn = NULL, Q = NULL, lmat = NULL,
                            psi.inv = NULL, c.data = NULL, ...) {
    load.psi    <- lapply(seq_along(nn), function(g) lmat[,,g] * psi.inv[,g])
    f.omega     <- lapply(seq_along(nn), function(g) diag(Q) + crossprod(load.psi[[g]], lmat[,,g]))
    U.f         <- lapply(f.omega, chol)
    z.f         <- lapply(seq_along(nn), function(g) matrix(rnorm(Q * nn[g], 0, 1), nr=Q, nc=nn[g]))
    v.f         <- lapply(seq_along(nn), function(g) backsolve(U.f[[g]], z.f[[g]]))
    mu.f        <- lapply(seq_along(nn), function(g) c.data[[g]] %*% (load.psi[[g]] %*% chol2inv(U.f[[g]])))
      lapply(seq_along(nn), function(g) mu.f[[g]] + t(v.f[[g]]))
  }

# Loadings
  sim.load      <- function(l.sigma = NULL, Q = NULL, c.data.j = NULL, f = NULL, 
                            psi.inv.j = NULL, FtF = NULL, G = NULL, ...) {
    l.omega     <- lapply(seq_len(G), function(g) l.sigma + psi.inv.j[g] * FtF[[g]])
    U.load      <- lapply(l.omega, chol)
    z.load      <- matrix(rnorm(Q * G, 0, 1), nr=Q, nc=G)
    v.load      <- do.call(cbind, lapply(seq_len(G), function(g) backsolve(U.load[[g]], z.load[,g])))
    mu.load     <- do.call(cbind, lapply(seq_len(G), function(g) psi.inv.j[g] * chol2inv(U.load[[g]]) %*% crossprod(f[[g]], c.data.j[[g]])))
      mu.load + v.load
  }

# Uniquenesses
  sim.psi.inv   <- function(N = NULL, P = NULL, psi.alpha = NULL, psi.beta = NULL, 
                            c.data = NULL, f = NULL, lmat = NULL, G = NULL, ...) { 
    rate.t      <- lapply(seq_len(G), function(g) c.data[[g]] - tcrossprod(f[[g]], lmat[,,g]))
    rate.t      <- unlist(lapply(seq_len(G), function(g) colSums(rate.t[[g]] * rate.t[[g]])))
      matrix(rgamma(P * G, shape=(N + psi.alpha)/2, 
                    rate=(rate.t + psi.beta)/2), nr=P, nc=G)
  }

# Priors
  # Means
    sim.mu.p    <- function(P = NULL, sigma.mu = NULL, G = NULL, ...) {
      U.mu      <- sqrt(1/sigma.mu)
      z.mu      <- matrix(rnorm(P * G, 0, 1), nr=P, ncol=G)
        U.mu * z.mu
    }
  
  # Scores
    sim.f.p     <- function(Q = NULL, nn = NULL, ...) {
        lapply(seq_along(nn), function(g) matrix(rnorm(nn[g] * Q, 0, 1), nr=nn[g], nc=Q))
    }

  # Loadings
    sim.load.p  <- function(Q = NULL, P = NULL, sigma.l = NULL, G = NULL, ...) {
      U.load    <- sqrt(1/sigma.l)
      z.load    <- array(rnorm(P * Q * G, 0, 1), dim=c(P, Q, G))
        z.load * U.load
    }

  # Uniquenesses
    sim.psi.p   <- function(P = NULL, psi.alpha = NULL, psi.beta = NULL, G = NULL, ...) {
        matrix(rgamma(n=P * G, shape=psi.alpha/2, rate=psi.beta/2), nr=P, nc=G) 
    }

  # Multivariate Normal Density
    mvdnorm     <- function(data = NULL, mu = NULL, Sigma = NULL, 
                            P = NULL, log.d = T, ...) {
      if(all(Sigma[!diag(P)] == 0)) {
        U.Sig   <- sqrt(Sigma)
      } else {
        U.Sig   <- chol(Sigma)
      }
      rooti     <- backsolve(U.Sig, diag(P))
      quad      <- crossprod(rooti, t(data) - mu)
      quads     <- colSums(quad * quad)
      log.dens  <- - P/2 * log(2 * pi) + sum(log(diag(rooti))) - 0.5 * quads
        if(isTRUE(log.d)) {
          log.dens
        } else {
          exp(log.dens)
        }
    }

  # Mixing Proportions
    sim.pi      <- function(alpha.pi = NULL, nn = 0, ...) {
        rdirichlet(1, alpha.pi + nn)
    }
  
  # Cluster Labels
    sim.z.p     <- function(N = NULL, prob.z = NULL, ...) {
      ind.mat   <- rmultinom(N, size=1, prob=prob.z)
      labs      <- which(ind.mat != 0, arr.ind=T)[,1]
        factor(labs, levels=seq_along(prob.z))
    }
