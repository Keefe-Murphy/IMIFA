###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals (Group Case)

  # Means
    sim.mu.m    <- function(nn = NULL, P = NULL, sigma.mu = NULL, psi.inv = NULL, 
                            sum.data = NULL, sum.f = NULL, lmat = NULL, G = NULL, ...) {
      mu.omega  <- 1/(sigma.mu + sweep(psi.inv, 2, nn, FUN="*"))
      U.mu      <- apply(mu.omega, 2, sqrt)
      z.mu      <- matrix(rnorm(P * G, 0, 1), nr=P, nc=G)
      v.mu      <- U.mu * z.mu
      lf.prod   <- do.call(cbind, lapply(seq_len(G), function(g) adrop(lmat[,,g, drop=F], drop=3) %*% sum.f[,g]))
      mu.mu     <- mu.omega * (psi.inv * (sum.data - lf.prod))
        mu.mu + v.mu
    }
  
  # Scores
    sim.score.m <- function(nn = NULL, Q = NULL, lmat = NULL,
                            psi.inv = NULL, c.data = NULL, ...) {
      load.psi  <- lapply(seq_along(nn), function(g) adrop(lmat[,,g, drop=F], drop=3) * psi.inv[,g])
      U.f       <- lapply(seq_along(nn), function(g) chol(diag(Q) + crossprod(load.psi[[g]], adrop(lmat[,,g, drop=F], drop=3))))
      z.f       <- lapply(seq_along(nn), function(g) matrix(rnorm(Q * nn[g], 0, 1), nr=Q, nc=nn[g]))
      v.f       <- lapply(seq_along(nn), function(g) backsolve(U.f[[g]], z.f[[g]]))
      mu.f      <- lapply(seq_along(nn), function(g) c.data[[g]] %*% (load.psi[[g]] %*% chol2inv(U.f[[g]])))
        do.call(rbind, lapply(seq_along(nn), function(g) mu.f[[g]] + t(v.f[[g]])))
    }
  
  # Loadings
    sim.load.m  <- function(l.sigma = NULL, Q = NULL, c.data.j = NULL, f = NULL, 
                            psi.inv.j = NULL, FtF = NULL, G = NULL, z = NULL, ...) {
      U.load    <- lapply(seq_len(G), function(g) chol(l.sigma + psi.inv.j[g] * FtF[[g]]))
      z.load    <- matrix(rnorm(Q * G, 0, 1), nr=Q, nc=G)
      v.load    <- do.call(cbind, lapply(seq_len(G), function(g) backsolve(U.load[[g]], z.load[,g])))
      mu.load   <- do.call(cbind, lapply(seq_len(G), function(g) psi.inv.j[g] * chol2inv(U.load[[g]]) %*% crossprod(f[z == g,, drop=F], c.data.j[[g]])))
        mu.load + v.load
    }
  
  # Uniquenesses
    sim.psi.im  <- function(N = NULL, P = NULL, psi.alpha = NULL, psi.beta = NULL, 
                            c.data = NULL, f = NULL, lmat = NULL, G = NULL, z = NULL, ...) { 
      rate.t    <- lapply(seq_len(G), function(g) c.data[[g]] - tcrossprod(f[z == g,, drop=F], adrop(lmat[,,g, drop=F], drop=3)))
      rate.t    <- unlist(lapply(rate.t, function(x) colSums(x * x)))
        matrix(rgamma(P * G, shape=(N + psi.alpha)/2, 
                      rate=(rate.t + psi.beta)/2), nr=P, nc=G)
    }
  
  # Mixing Proportions
    sim.pi      <- function(pi.alpha = NULL, nn = 0, ...) {
        rdirichlet(1, pi.alpha + nn)
    }
    
  # Cluster Labels
    sim.z       <- function(data = NULL, mu = NULL, Sigma = NULL, N = NULL, 
                            G = NULL, P = NULL, pi.prop = NULL, ...) {
      numer     <- do.call(cbind, lapply(seq_len(G), function(g) exp(mvdnorm(data, mu[,g], Sigma[[g]], P, log.d=T) + log(pi.prop[,g]))))
      denomin   <- rowSums(numer)
      pz        <- sweep(numer, 1, denomin, FUN="/")
      pz[rowSums(pz > 0) == 0,] <- rep(1/G, G)
      pz[pz <= 0]               <- .Machine$double.eps
      z         <- factor(do.call(c, lapply(seq_len(N), function(i) which(rmultinom(1, size=1, prob=pz[i,]) !=0))), levels=seq_len(G))
        return(list(z = z, log.liks = log(denomin)))
    }

# Priors (Group Case)

  # Means
    sim.mu.mp   <- function(P = NULL, sigma.mu = NULL, G = NULL, ...) {
      U.mu      <- sqrt(1/sigma.mu)
      z.mu      <- matrix(rnorm(P * G, 0, 1), nr=P, nc=G)
        U.mu * z.mu
    }
  
  # Scores
    sim.f.mp    <- function(Q = NULL, nn = NULL, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }

  # Loadings
    sim.load.mp <- function(Q = NULL, P = NULL, sigma.l = NULL, G = NULL, ...) {
      U.load    <- sqrt(1/sigma.l)
      z.load    <- array(rnorm(P * Q * G, 0, 1), dim=c(P, Q, G))
        z.load * U.load
    }

  # Uniquenesses
    sim.psi.imp <- function(P = NULL, psi.alpha = NULL, psi.beta = NULL, G = NULL, ...) {
        matrix(rgamma(n=P * G, shape=psi.alpha/2, rate=psi.beta/2), nr=P, nc=G) 
    }

  # Cluster Labels
    sim.z.p     <- function(N = NULL, prob.z = NULL, ...) {
      ind.mat   <- rmultinom(N, size=1, prob=prob.z)
      labs      <- which(ind.mat != 0, arr.ind=T)[,1]
        factor(labs, levels=seq_along(prob.z))
    }

# Full Conditionals (Single Case)

  # Means
    sim.mu      <- function(N = NULL, P = NULL, sigma.mu = NULL, psi.inv = NULL, 
                              sum.data = NULL, sum.f = NULL, lmat = NULL, ...) {
      mu.omega  <- 1/(sigma.mu + N * psi.inv)
      U.mu      <- sqrt(mu.omega)
      z.mu      <- rnorm(P, 0, 1)
      v.mu      <- U.mu * z.mu
      mu.mu     <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f))
        as.vector(mu.mu + v.mu)
    }
  
  # Scores
    sim.score   <- function(N = NULL, Q = NULL, lmat = NULL,
                              psi.inv = NULL, c.data = NULL, ...) {
      load.psi  <- lmat * psi.inv
      U.f       <- chol(diag(Q) + crossprod(load.psi, lmat))
      z.f       <- matrix(rnorm(Q * N, 0, 1), nr=Q, nc=N)
      v.f       <- backsolve(U.f, z.f)
      mu.f      <- c.data %*% (load.psi %*% chol2inv(U.f))
        mu.f + t(v.f)
    }
  
  # Loadings
    sim.load    <- function(l.sigma = NULL, Q = NULL, c.data.j = NULL, 
                            f = NULL, psi.inv.j = NULL, FtF = NULL, ...) {
      U.load    <- chol(l.sigma + psi.inv.j * FtF)
      z.load    <- rnorm(Q, 0, 1)
      v.load    <- backsolve(U.load, z.load)
      mu.load   <- psi.inv.j * chol2inv(U.load) %*% crossprod(f, c.data.j)
        mu.load + v.load
    }
  
  # Uniquenesses
    sim.psi.i   <- function(N = NULL, P = NULL, psi.alpha = NULL, psi.beta = NULL, 
                            c.data = NULL, f = NULL, lmat = NULL, ...) { 
      rate.t    <- c.data - tcrossprod(f, lmat)
      rate.t    <- colSums(rate.t * rate.t)
        rgamma(P, shape=(N + psi.alpha)/2, 
               rate=(rate.t + psi.beta)/2) 
    }

  # Local Shrinkage
    sim.phi     <- function(Q = NULL, P = NULL, phi.nu = NULL, 
                            tau = NULL, load.2 = NULL, ...) {
      rate.t    <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
        matrix(rgamma(P * Q, shape=(phi.nu + 1)/2, 
                      rate=rate.t), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta1  <- function(Q = NULL, P = NULL, alpha.d1 = NULL, delta = NULL,
                            tau = NULL, sum.term = NULL, ...) {
        rgamma(1, shape=alpha.d1 + P * Q/2, 
               rate=1 + 0.5/delta[1] * tau %*% sum.term)
    }
    
    sim.deltak  <- function(Q = NULL, P = NULL, k = NULL, alpha.d2 = NULL,
                            delta = NULL, tau = NULL, sum.term = NULL) {
        rgamma(1, shape=alpha.d2 + P/2 * (Q - k + 1), 
               rate=1 + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
    }

# Priors (Single Case)

  # Means
    sim.mu.p    <- function(P = NULL, sigma.mu = NULL, ...) {
      U.mu      <- sqrt(1/sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu
    }
  
  # Scores
    sim.f.p     <- function(Q = NULL, N = NULL, ...) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }
  
  # Loadings
    sim.load.p  <- function(Q = NULL, P = NULL, sigma.l = NULL, ...) {
      U.load    <- sqrt(1/sigma.l)
      z.load    <- matrix(rnorm(P * Q, 0, 1), nr=P, nc=Q)
        z.load * U.load
    }

  # Loadings (Shrinkage)
    sim.load.ps <- function(Q = NULL, D.load = NULL, ...) {
      U.load    <- sqrt(1/D.load)
      z.load    <- rnorm(Q, 0, 1)
        U.load * z.load
    }
  
  # Uniquenesses
    sim.psi.ip  <- function(P = NULL, psi.alpha = NULL, psi.beta = NULL, ...) {
        rgamma(n=P, shape=psi.alpha/2, rate=psi.beta/2) 
    }

  # Local Shrinkage
    sim.phi.p   <- function(Q = NULL, P = NULL, phi.nu = NULL, ...) {
        matrix(rgamma(n=P * Q, shape=phi.nu/2, rate=phi.nu/2), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta.p <- function(Q = NULL, alpha.d1 = NULL, alpha.d2 = NULL, ...) {
      delta1    <- rgamma(n=1,     shape=alpha.d1, rate=1)
      deltak    <- rgamma(n=Q - 1, shape=alpha.d2, rate=1)
        c(delta1, deltak)
    }

# Other Functions

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