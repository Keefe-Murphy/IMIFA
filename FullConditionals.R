###############################
### IMIFA Full Conditionals ###
###############################

# Full Conditionals

  # Means
    sim.mu      <- function(N, P, mu.sigma, psi.inv, sum.data, sum.f, lmat, mu.zero) {
      mu.omega  <- 1/(mu.sigma + N * psi.inv)
        as.vector(mu.omega * (psi.inv * (sum.data - lmat %*% sum.f) + mu.sigma * mu.zero) + sqrt(mu.omega) * rnorm(P))
    }
  
  # Scores
    sim.score   <- function(N, Q, lmat, psi.inv, c.data, Q1) {
      load.psi  <- lmat * psi.inv
      U.f       <- diag(Q) + crossprod(load.psi, lmat)
      U.f       <- if(Q1) chol(U.f) else sqrt(U.f)
      mu.f      <- c.data %*% (load.psi %*% if(Q1) chol2inv(U.f) else 1/(U.f * U.f))
        mu.f + t(backsolve(U.f, matrix(rnorm(Q * N), nr=Q, nc=N)))
    }
      
  # Loadings
    sim.load    <- function(l.sigma, Q, c.data, P, f, phi, tau, psi.inv, FtF, Q1, shrink = T) {
      U.load    <- if(shrink) phi * tau * diag(Q) + psi.inv * FtF else l.sigma + psi.inv * FtF
      U.load    <- if(Q1) chol(U.load) else sqrt(U.load)
      mu.load   <- psi.inv * (if(Q1) chol2inv(U.load) else 1/(U.load * U.load)) %*% crossprod(f, c.data)
        as.vector(mu.load + backsolve(U.load, rnorm(Q)))
    }
    
  # Uniquenesses
    sim.psi.i   <- function(N, P, psi.alpha, psi.beta, c.data, f, lmat) { 
      rate.t    <- c.data - tcrossprod(f, lmat)
        rgamma(P, shape=N/2 + psi.alpha, rate=colSums(rate.t * rate.t)/2 + psi.beta) 
    }

  # Local Shrinkage
    sim.phi     <- function(Q, P, phi.nu, tau, load.2) {
        matrix(rgamma(P * Q, shape=1/2 + phi.nu, rate=(phi.nu + sweep(load.2, 2, tau, FUN="*"))/2), nr=P, nc=Q)
    }
  
  # Global Shrinkage
    sim.delta1  <- function(Q, P, alpha.d1, delta.1, beta.d1, tau, sum.term) {
        rgamma(1, shape=alpha.d1 + P * Q/2, rate=beta.d1 + 0.5/delta.1 * tau %*% sum.term)
    }
    
    sim.deltak  <- function(Q, P, k, alpha.dk, beta.dk, delta.k, tau.kq, sum.term.kq) {
        rgamma(1, shape=alpha.dk + P/2 * (Q - k + 1), rate=beta.dk + 0.5/delta.k * tau.kq %*% sum.term.kq)
    }

  # Mixing Proportions
    sim.pi      <- function(pi.alpha, nn, inf.G=F, len) {
      if(inf.G) {
        vs      <- rbeta(length(nn), 1 + nn, pi.alpha + sum(nn) - cumsum(nn))
        vs[len] <- 1
          return(list(Vs = vs, pi.prop = do.call(cbind, lapply(seq_along(nn), function(t) vs[t] * prod((1 - vs[seq_len(t - 1)]))))))
      } else {
          rdirichlet(1, pi.alpha + nn)
      }
    }
  
  # Cluster Labels
    sim.z       <- function(data, mu, Sigma, Gseq, N, pi.prop, slice.ind=NULL, Q0) {
      log.num   <- do.call(cbind, lapply(Gseq, function(g, q=Q0[g]) dmvn(data, mu[,g], if(q) Sigma[[g]] else sqrt(Sigma[[g]]), log=T, isChol=!q) + log(pi.prop[,g])))
      if(!missing(slice.ind)) {
        log.num <- log.num + log(slice.ind)
      }
      log.denom <- rowLogSumExps(log.num)
      lnp       <- sweep(log.num, 1, log.denom, FUN="-")
      for(g in Gseq[-1])      {
        lnp[,g] <- colLogSumExps(rbind(lnp[,g], lnp[,g - 1]))
      }
        return(list(z = unname(rowSums(-rexp(N) > lnp) + 1), log.likes = log.denom))
    }

  # Alpha
    sim.alpha   <- function(lower, upper, trunc.G, alpha, Vs) {
      alpha.new <- runif(1, lower, upper)
      a.prob    <- trunc.G * (log(alpha.new) - log(alpha)) + (alpha.new - alpha) * sum(log((1 - Vs[-trunc.G])))
      accept    <- a.prob >= 0 || -rexp(1) < a.prob
        return(list(alpha = ifelse(accept, alpha.new, alpha), rate = accept))
    }

# Priors

  # Means
    sim.mu.p    <- function(P, sigma.mu, mu.zero) {
      sqrt(sigma.mu) * rnorm(P) + mu.zero
    }
  
  # Scores
    sim.f.p     <- function(Q, N) {
        matrix(rnorm(N * Q), nr=N, nc=Q)
    }
  
  # Loadings
    sim.load.p  <- function(Q, P, sigma.l, phi, tau, shrink = T) {
      if(shrink) {
        sqrt(1/(phi * tau)) * rnorm(Q)
      } else     {
        sqrt(sigma.l) * matrix(rnorm(P * Q), nr=P, nc=Q)
      }
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
        factor(which(rmultinom(N, size=1, prob=prob.z) != 0, arr.ind=T)[,1], levels=seq_along(prob.z))
    }

# Other Functions
    
  # Uniqueness Hyperparameters
    psi.hyper   <- function(alpha, covar) {
      inv.cov   <- try(chol2inv(chol(covar)), silent=T)
      if(inherits(inv.cov, "try-error"))   {
        inv.cov <- 1/covar
      }
        (alpha - 1)/diag(inv.cov) 
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
        z.perm  <- c(z.perm, setNames(setdiff(Gs, z.perm), setdiff(Gs, names(z.perm))))
      }
      z.names   <- as.numeric(names(z.perm))
      z.perm    <- z.perm[order(z.names)]
      z.sw      <- factor(z.new, labels=z.perm[seq_along(ng[ng > 0])])
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
          if(all(len == range.G)) obj0g <- if(switch0g) lapply(seq_along(range.G), function(g) matrix(obj0g[[g]], nr=1)) else stop(paste0(sw.name, " must be TRUE if the dimension of ", obj.name, " depends on G"))
          if(all(len == V))       obj0g <- lapply(seq_along(range.G), function(g) matrix(obj0g[[g]], nr=V, nc=range.G[g]))
        } else if(!all(sapply(seq_along(range.G), function(g) is.matrix(obj0g[[g]]) && any(identical(dim(obj0g[[g]]), c(1, range.G[g])), identical(dim(obj0g[[g]]), c(V, range.G[g])))))) {
                                           stop(paste0(ifelse(length(range.G) > 1, "Each element of ", ""), obj.name, " must be either of length 1, ", ifelse(P.dim, paste0("P=", V, ", or it's corresponding range.G, or a matrix with P rows and it's corresponding range.G columns"), paste0("or G=", range.G)))) 
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
                            byrank=F, breaks=length(colors)) { 
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