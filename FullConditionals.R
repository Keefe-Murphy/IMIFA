###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    sim.mu      <- function(N, P, mu.sigma, psi.inv, sum.data, sum.f, lmat, mu.zero) {
      mu.omega  <- 1/(mu.sigma + N * psi.inv)
      U.mu      <- sqrt(mu.omega)
      z.mu      <- rnorm(P, 0, 1)
      v.mu      <- U.mu * z.mu
      mu.mu     <- mu.omega * (psi.inv * (sum.data - lmat %*% sum.f) + mu.sigma * mu.zero)
        as.vector(mu.mu + v.mu)
    }
  
  # Scores
    sim.score   <- function(N, Q, lmat, psi.inv, c.data) {
      load.psi  <- lmat * psi.inv
      U.f       <- chol(diag(Q) + crossprod(load.psi, lmat))
      z.f       <- matrix(rnorm(Q * N, 0, 1), nr=Q, nc=N)
      v.f       <- backsolve(U.f, z.f)
      mu.f      <- c.data %*% (load.psi %*% chol2inv(U.f))
        mu.f + t(v.f)
    }
  
  # Loadings
    sim.load    <- function(l.sigma, Q, c.data, P, f, phi, tau, psi.inv, FtF, shrink = T) {
      if(shrink) {
        U.load  <- chol(phi * tau * diag(Q) + psi.inv * FtF)
      } else     {
        U.load  <- chol(l.sigma + psi.inv * FtF)
      }
      z.load    <- rnorm(Q, 0, 1)
      v.load    <- backsolve(U.load, z.load)
      mu.load   <- psi.inv * chol2inv(U.load) %*% crossprod(f, c.data)
        as.vector(mu.load + v.load)
    }

  # Uniquenesses
    sim.psi.i   <- function(N, P, psi.alpha, psi.beta, c.data, f, lmat) { 
      rate.t    <- c.data - tcrossprod(f, lmat)
      rate.t    <- colSums(rate.t * rate.t)
        rgamma(P, shape=N/2 + psi.alpha, 
               rate=rate.t/2 + psi.beta) 
    }

  # Local Shrinkage
    sim.phi     <- function(Q, P, phi.nu, tau, load.2) {
      rate.t    <- (phi.nu + sweep(load.2, 2, tau, FUN="*"))/2
        matrix(rgamma(P * Q, shape=1/2 + phi.nu, 
                      rate=rate.t), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta1  <- function(Q, P, alpha.d1, delta, beta.d1, tau, sum.term) {
        rgamma(1, shape=alpha.d1 + P * Q/2, 
               rate=beta.d1 + 0.5/delta[1] * tau %*% sum.term)
    }
    
    sim.deltak  <- function(Q, P, k, alpha.dk, beta.dk, delta, tau, sum.term) {
        rgamma(1, shape=alpha.dk + P/2 * (Q - k + 1), 
               rate=beta.dk + 0.5/delta[k] * tau[k:Q] %*% sum.term[k:Q])
    }

  # Mixing Proportions
    sim.pi      <- function(pi.alpha, nn, inf.G=F) {
      if(inf.G) {
        vs      <- rbeta(length(nn), 1 + nn, pi.alpha - cumsum(nn))
        vs[length(vs)]    <- 1
          do.call(cbind, lapply(seq_along(nn), function(t) vs[t] * prod((1 - vs[seq_len(t - 1)]))))  
      } else {
          rdirichlet(1, pi.alpha + nn)
      }
    }
  
  # Cluster Labels
    sim.z       <- function(data, mu, Sigma, Gseq, N, pi.prop, slice.ind=NULL) {
      log.num   <- do.call(cbind, lapply(Gseq, function(g) mvdnorm(data, mu[,g], Sigma[[g]], log.d=T) + log(pi.prop[,g])))
      if(!missing(slice.ind)) {
        log.num <- log.num + log(slice.ind)
      }
      log.denom <- rowLogSumExps(log.num)
      lnp       <- sweep(log.num, 1, log.denom, FUN="-")
      for(g in Gseq[-1])      {
        lnp[,g] <- colLogSumExps(rbind(lnp[,g], lnp[,g - 1]))
      }
      z         <- rowSums(-rexp(N) > lnp) + 1
        return(list(z = unname(z), log.likes = log.denom))
    }

# Priors

  # Means
    sim.mu.p    <- function(P, sigma.mu, mu.zero) {
      U.mu      <- sqrt(sigma.mu)
      z.mu      <- rnorm(P, 0, 1)
        U.mu * z.mu + mu.zero
    }
  
  # Scores
    sim.f.p     <- function(Q, N) {
        matrix(rnorm(N * Q, 0, 1), nr=N, nc=Q)
    }
  
  # Loadings
    sim.load.p  <- function(Q, P, sigma.l, phi, tau, shrink = T) {
      if(shrink) {
        U.load  <- sqrt(1/(phi * tau))
        z.load  <- rnorm(Q, 0, 1)
      } else     {
        U.load  <- sqrt(sigma.l)
        z.load  <- matrix(rnorm(P * Q, 0, 1), nr=P, nc=Q)
      }
        U.load * z.load
    }
  
  # Uniquenesses
    sim.psi.ip  <- function(P, psi.alpha, psi.beta) {
        rgamma(n=P, shape=psi.alpha, rate=psi.beta) 
    }

  # Local Shrinkage
    sim.phi.p   <- function(Q, P, phi.nu) {
        matrix(rgamma(n=P * Q, shape=phi.nu, rate=phi.nu), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta.p <- function(Q=2, alpha, beta) {
        rgamma(n=Q - 1, shape=alpha, rate=beta)
    }
        
  # Cluster Labels
    sim.z.p     <- function(N, prob.z) {
      ind.mat   <- rmultinom(N, size=1, prob=prob.z)
      labs      <- which(ind.mat != 0, arr.ind=T)[,1]
        factor(labs, levels=seq_along(prob.z))
    }

# Other Functions

  # Multivariate Normal Density
    mvdnorm     <- function(data, mu, Sigma, log.d = T) {
      P         <- length(mu)
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
  
  # Uniqueness Hyperparameters
    psi.hyper   <- function(alpha, covar) {
      inv.cov   <- try(chol2inv(chol(covar)), silent=T)
      if(inherits(inv.cov, "try-error"))   {
        inv.cov <- 1/covar
      }
      psi.beta  <- (alpha - 1)/diag(inv.cov) 
        psi.beta
    }

  # Label Switching
    lab.switch  <- function(z.new, z.old, Gs, ng=tabulate(z.new)) {
      tab       <- table(z.new, z.old, dnn=NULL)
      tab.tmp   <- tab[rowSums(tab) != 0,colSums(tab) != 0, drop=F]
      nc        <- ncol(tab.tmp)
      nr        <- nrow(tab.tmp)
      if(nc > nr) {
        tmp.mat <- matrix(rep(0, nc), nr=nc - nr, nc=nc)
        rownames(tmp.mat) <- setdiff(as.numeric(colnames(tab.tmp)), as.numeric(rownames(tab.tmp)))[seq_len(nc - nr)]
        tab.tmp <- rbind(tab.tmp, tmp.mat)
        tab.tmp <- tab.tmp[match(rownames(tab.tmp), Gs),]
      } else if(nr > nc) {
        tmp.mat <- matrix(rep(0, nr), nr=nr, nc=nr - nc)
        colnames(tmp.mat) <- setdiff(as.numeric(rownames(tab.tmp)), as.numeric(colnames(tab.tmp)))[seq_len(nr - nc)]
        tab.tmp <- cbind(tab.tmp, tmp.mat)
        tab.tmp <- tab.tmp[,match(colnames(tab.tmp), Gs)]
      }
      if(nr == 1) {
        z.perm  <- setNames(as.numeric(colnames(tab.tmp)), as.numeric(colnames(tab.tmp)))
      } else if(nc == 1) {
        z.perm  <- setNames(as.numeric(colnames(tab.tmp)), as.numeric(colnames(tab.tmp)))
      } else {
        z.perm  <- suppressWarnings(matchClasses(tab.tmp, method="exact", verbose=F))
        z.perm  <- setNames(as.numeric(z.perm), names(z.perm))
      }
      if(length(Gs) > length(z.perm)) {
        miss.z1 <- setdiff(Gs, names(z.perm))
        miss.z2 <- setdiff(Gs, z.perm)
        z.perm  <- c(z.perm, setNames(miss.z2, miss.z1))
      }
      z.names   <- as.numeric(names(z.perm))
      z.perm    <- z.perm[order(z.names)]
      ng        <- ng[ng > 0]
      z.sw      <- factor(z.new, labels=z.perm[seq_along(ng)])
        return(list(z = as.numeric(levels(z.sw))[z.sw], z.perm = z.perm))
    }

  # Length Checker
    len.check   <- function(obj0g, switch0g, P.dim=T) {
      V         <- ifelse(P.dim, P, 1)
      obj.name  <- deparse(substitute(obj0g))
      sw.name   <- deparse(substitute(switch0g))
      if(!is.list(obj0g))        obj0g  <- list(obj0g)
      if(length(obj0g) != length(range.G)) stop(paste0(obj.name, " must be a list of length ", length(range.G)))
      len       <- sapply(obj0g, length)
      if(is.element(method, c("FA", "IFA"))) {
        if(any(!is.element(len, c(1, V)))) stop(paste0(obj.name, " must be list of length 1 containing a scalar", ifelse(P.dim, paste0(" or a vector of length P=", V), ""), " for a 1-group model"))
      } else {
        if(all(is.element(len, c(1, range.G, V)))) {
          if(all(len == 1))       obj0g <- lapply(seq_along(range.G), function(g) matrix(obj0g[[g]], nr=1, nc=range.G[g]))
          if(all(len == range.G)) obj0g <- if(switch0g) lapply(seq_along(range.G), function(g) matrix(obj0g[[g]], nr=1)) else stop(paste0(sw.name, "must be TRUE if the dimension of ", obj.name, " depends on G"))
          if(all(len == V))       obj0g <- lapply(seq_along(range.G), function(g) matrix(obj0g[[g]], nr=V, nc=range.G[g]))
        } else if(!all(sapply(seq_along(range.G), function(g) is.matrix(obj0g[[g]]) && any(identical(dim(obj0g[[g]]), c(1, range.G[g])), identical(dim(obj0g[[g]]), c(V, range.G[g])))))) {
                                           stop(paste0("Each element of ", obj.name, " must be either of length 1, P=", V, ", or it's corresponding range.G, or a matrix with P rows and it's corresponding range.G columns")) 
        } else if(all(sapply(obj0g, is.matrix), !switch0g) && any(sapply(seq_along(range.G), function(g) any(dim(obj0g[[g]]) == range.G[g])))) {
                                           stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
        }
      }
      if(all(length(unique(unlist(obj0g))) > 1,
             !switch0g, !P.dim))           stop(paste0(obj.name, " must be a scalar if ", sw.name, " is TRUE"))
        obj0g
    }

  # Loadings Heatmaps
    mat2color   <- function(m, colors = heat.colors(30), 
                            byrank=FALSE, breaks=length(colors)) { 
      m1        <- if(byrank == T) rank(m) else m
      facs      <- cut(m1, breaks, include.lowest=TRUE)
      answer    <- colors[as.numeric(facs)]
      if(is.matrix(m)) {
        answer  <- matrix(answer, nrow(m), ncol(m))
        rownames(answer)  <- rownames(m)
        colnames(answer)  <- colnames(m)
      }
        answer    
    }