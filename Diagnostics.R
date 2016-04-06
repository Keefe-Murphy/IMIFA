#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.imifa       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL, Q = NULL, Q.meth = c("Mode", "Median"),
                             criterion = c("bicm", "aicm", "bic.mcmc", "aic.mcmc"), recomp = F, ...) {
  defpar         <- par(no.readonly = T)
  defop          <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(sims))               stop("Simulations must be supplied")
  if(!exists(deparse(substitute(sims)),
             envir=.GlobalEnv))   stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")      stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  store          <- seq(from=burnin + 1, to=attr(sims, "Store"), by=thinning)
  n.store        <- length(store)
  if(n.store     <= 1)            stop(paste0("burnin must be less than the stored number of iterations"))
  method         <- attr(sims, "Method")
  n.fac          <- attr(sims, "Factors")
  n.grp          <- attr(sims, "Groups")
  n.obs          <- attr(sims, "Obs")
  n.var          <- attr(sims, "Vars")
  sw             <- attr(sims, "Switch")
  criterion      <- match.arg(criterion)
  cent           <- attr(sims, "Center")
  scaling        <- attr(sims, "Scaling")
  if(all(method  == "MIFA", 
     !is.element(criterion, 
     c("aicm", "bicm"))))         stop("criterion should be one of 'aicm' or 'bicm' for the MIFA method")
  if(!is.logical(recomp))         stop("recomp must be TRUE or FALSE")
  if(any(burnin   > 0, 
     thinning     > 1)) recomp <- T
  
  G.ind          <- 1
  Q.ind          <- 1
  G.T            <- !missing(G)
  Q.T            <- !missing(Q)
  if(G.T) {
    if(!is.element(method, c("FA", "IFA"))) {
      G.ind      <- which(n.grp == G)
    } else if(G > 1)              warning(paste0("G must be equal to 1 for the ", method, " method"), call.=F)
    if(all(is.element(method, c("MFA", "MIFA")),
       !is.element(G, n.grp)))    stop("This G value was not used during simulation")
  } 
  if(Q.T) {
    if(!is.element(method, c("IFA", "MIFA", "IMIFA"))) {
      Q.ind      <- which(n.fac == Q)
    }
    if(all(is.element(method, c("FA", "MFA")),
       !is.element(Q, n.fac)))    stop("This Q value was not used during simulation")
    if(all(is.element(method, c("IFA", "classify")), 
      (Q * (n.fac - Q)) < 0))     stop(paste0("Q cannot be greater than the number of factors in ", match.call()$sims))
  } 
  G              <- ifelse(all(G.T, !is.element(method, c("FA", "IFA"))), G, 1)
  
  if(method == "IFA") {
    if(missing(Q.meth)) {
      Q.meth     <- "Mode"
    } else   {
      Q.meth     <- match.arg(Q.meth)
    }
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.store      <- sims[[G.ind]][[Q.ind]]$Q.store[store]
    Q.tab        <- table(Q.store, dnn=NULL)
    Q.prob       <- prop.table(Q.tab)
    Q.mode       <- as.numeric(names(Q.tab[Q.tab == max(Q.tab)]))
    Q.med        <- ceiling(median(Q.store) * 2)/2
    if(!Q.T) {
      if(Q.meth  == 'Mode') { 
        Q        <- min(Q.mode)
      } else {
        Q        <- Q.med
      }
    }
    Q.CI         <- round(quantile(Q.store, c(0.025, 0.975)))
    GQ.res       <- list(G = G, Q = Q, Mode = Q.mode, Median = Q.med, 
                         CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
  }
    
  if(is.element(method, c("FA", "MFA"))) {
    G.range      <- length(n.grp)
    Q.range      <- length(n.fac)
    
  # Retrieve log-likelihoods and tune G & Q according to criterion
    aicm         <- bicm       <- 
    aic.mcmc     <- bic.mcmc   <- matrix(NA, nr=G.range, nc=Q.range, dimnames=list(paste0("G", n.grp), paste0("Q", n.fac)))
    for(g in seq_len(G.range)) { 
      for(q in seq_len(Q.range)) {
       log.likes     <- sims[[g]][[q]]$ll.store[store]
       K             <- attr(sims[[g]][[q]], "K")
       ll.max        <- 2 * max(log.likes)
       ll.var        <- 2 * var(log.likes)
       ll.mean       <- mean(log.likes)
       aicm[g,q]     <- ll.max - ll.var * 2
       bicm[g,q]     <- ll.max - ll.var * log(n.obs)     
       aic.mcmc[g,q] <- ll.max - K * 2
       bic.mcmc[g,q] <- ll.max - K * log(n.obs)
      }  
    }
    crit         <- get(criterion)
    crit.max     <- which(crit == max(crit), arr.ind = T)
  
  # Control for supplied values of G &/or Q
    if(all(Q.T, G.T)) {
      aic.mcmc   <- aic.mcmc[G.ind,Q.ind, drop=F]
      bic.mcmc   <- bic.mcmc[G.ind,Q.ind, drop=F]
      aicm       <- aicm[G.ind,Q.ind, drop=F]
      bicm       <- bicm[G.ind,Q.ind, drop=F]
    } else if(Q.T) {
      aic.mcmc   <- aic.mcmc[,Q.ind, drop=F]
      bic.mcmc   <- bic.mcmc[,Q.ind, drop=F]
      aicm       <- aicm[,Q.ind, drop=F]
      bicm       <- bicm[,Q.ind, drop=F]
      crit       <- crit[,Q.ind]
      G.ind      <- which(crit == max(crit))
      G          <- n.grp[G.ind]
    } else if(G.T) {
      aic.mcmc   <- aic.mcmc[G.ind,, drop=F]
      bic.mcmc   <- bic.mcmc[G.ind,, drop=F]
      aicm       <- aicm[G.ind,, drop=F]
      bicm       <- bicm[G.ind,, drop=F]
      crit       <- crit[G.ind,]
      Q.ind      <- which(crit == max(crit))
      Q          <- n.fac[Q.ind]
    } else {
      G.ind      <- crit.max[1]
      Q.ind      <- crit.max[2]
      G          <- n.grp[G.ind]
      Q          <- n.fac[Q.ind]
    }
    Q            <- setNames(rep(Q, G), paste0("Qg", seq_len(G)))
    GQ.res       <- list(G = G, Q = Q, AICM = aicm, BICM = bicm,
                         AIC.mcmc = aic.mcmc, BIC.mcmc = bic.mcmc)
  }
  
# Retrieve cluster labels and mixing proportions
  if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
    z            <- sims[[G.ind]][[Q.ind]]$z[,store]
    post.z       <- setNames(apply(z, 1, function(x) factor(which.max(tabulate(x)), levels=seq_len(G))), seq_len(n.obs))
    if(sw["pi.sw"])    {
      pi.prop    <- sims[[G.ind]][[Q.ind]]$pi.prop[,store]
      post.pi    <- rowMeans(pi.prop, dims=1)
    } else {
      post.pi    <- setNames(prop.table(tabulate(post.z, nbins=G)), paste0("Group ", seq_len(G)))
    }
    cluster      <- list(post.z = post.z, post.pi = post.pi, z = z)
    cluster      <- c(cluster, if(sw["pi.sw"]) list(pi.prop  = pi.prop))
    attr(cluster, "Z.init")    <- attr(sim[[G.ind]][[Q.ind]], "Z.init")
  }
  
# Retrieve (unrotated) scores
  if(all(Q == 0)) {
    if(sw["f.sw"])                warning("Scores not stored as model has zero factors", call.=F)
    sw["f.sw"]   <- F
  }
  if(sw["f.sw"])  {
    Qms          <- seq_len(max(Q))
    if(is.element(method, c("FA", "MFA"))) {
      f          <- sims[[G.ind]][[Q.ind]]$f[,Qms,store, drop=F]
    }
    if(is.element(method, c("IFA", "MIFA", "IMIFA"))) {
      f          <- as.array(sims[[G.ind]][[Q.ind]]$f)[,Qms,store, drop=F]
    }
  }

# Loop over g in G to extract other results
  result         <- list(list())
  temp.b         <- max(1, burnin)
  MSE  <- RMSE   <- NRMSE   <- CVRMSE   <- MAD   <- rep(NA, G)
  for(g in seq_len(G)) {
    Qg           <- Q[g]
    Qgs          <- seq_len(Qg)
    sw["l.sw"]   <- attr(sims, "Switch")["l.sw"]
    if(Qg == 0) {
      if(sw["l.sw"])              warning(paste0("Loadings not stored as", ifelse(G > 1, paste0(" group ", g), " model"), " has zero factors"), call.=F)
      sw["l.sw"] <- F
    }
  
  # Retrieve (unrotated) loadings  
    if(sw["l.sw"]) {
      if(all(method == "MFA", G > 1)) {
        lmat     <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,g,store, drop=F], drop=3)
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,g,temp.b, drop=F], drop=3:4)
      }
      if(any(method == "FA",  all(method == "MFA",  G == 1))) {
        lmat     <- sims[[G.ind]][[Q.ind]]$load[,Qgs,store, drop=F]
        l.temp   <- adrop(sims[[G.ind]][[Q.ind]]$load[,Qgs,temp.b, drop=F], drop=3)
      }
    }
  
    if(any(method == "IFA", all(method == "MIFA", G == 1))) {
      store      <- store[which(Q.store >= Qg)]
      n.store    <- length(store)
      temp.b     <- store[1]
      if(sw["l.sw"]) {
        lmat     <- as.array(sims[[G.ind]][[Q.ind]]$load)[,Qgs,store, drop=F]
        l.temp   <- adrop(as.array(sims[[G.ind]][[Q.ind]]$load)[,Qgs,temp.b, drop=F], drop=3)
      }
    }
    
  # Loadings matrix / identifiability / error metrics / etc.  
    if(sw["l.sw"])     {
      for(p in seq_len(n.store)) {
        rot          <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
        lmat[,,p]    <- adrop(lmat[,,p, drop=F], drop=3) %*% rot
        if(sw["f.sw"]) {
          if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
            f[post.z == g,,p]  <- adrop(f[post.z == g,,p, drop=F], drop=3) %*% rot
          } else {
            f[,,p]   <- adrop(f[,,p, drop=F], drop=3)    %*% rot
          }
          scores     <- list(f = f, post.f = rowMeans(f, dims=2))
        }  
      }
    }
  
  # Retrieve means, uniquenesses & empirical covariance matrix
    if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
      post.mu    <- sims[[G.ind]][[Q.ind]]$post.mu[,g]
      post.psi   <- sims[[G.ind]][[Q.ind]]$post.psi[,g]
      if(sw["mu.sw"])  {
        mu       <- sims[[G.ind]][[Q.ind]]$mu[,g,store]                            
      }
      if(sw["psi.sw"]) {
        psi      <- sims[[G.ind]][[Q.ind]]$psi[,g,store]
      }
      data       <- attr(sims, "Name")
      if(!exists(data,
         envir=.GlobalEnv)) {     warning(paste0("Object ", data, " not found in .GlobalEnv: can't compute empirical covariance and error metrics"), call.=F)
      } else {
        data     <- as.data.frame(get(data))
        data     <- data[sapply(data, is.numeric)]
        data     <- scale(data, center=cent, scale=scaling)
      }
      varnames   <- colnames(data)
      cov.emp    <- cov(data[post.z == g,, drop=F])
      dimnames(cov.emp)        <- list(varnames, varnames)
    } else {
      post.mu    <- sims[[G.ind]][[Q.ind]]$post.mu
      post.psi   <- sims[[G.ind]][[Q.ind]]$post.psi
      if(sw["mu.sw"])  {
        mu       <- sims[[G.ind]][[Q.ind]]$mu[,store]                            
      }
      if(sw["psi.sw"]) {
        psi      <- sims[[G.ind]][[Q.ind]]$psi[,store]
      }
      cov.emp    <- sims[[G.ind]][[Q.ind]]$cov.emp
    }
  
  # Compute posterior means and % variation explained
    if(sw["mu.sw"])  post.mu   <- rowMeans(mu, dims=1)
    if(sw["psi.sw"]) post.psi  <- rowMeans(psi, dims=1)
    if(sw["l.sw"]) { post.load <- rowMeans(lmat, dims=2)
              class(post.load) <- "loadings"
      var.exp    <- sum(colSums(post.load * post.load))/n.var
    } else   {
      var.exp    <- (sum(diag(cov.emp)) - sum(post.psi))/n.var
    }
  
  # Calculate estimated covariance matrices & compute error metrics
    if(all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)) {
      if(all(sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        if(Qg > 0) {
          cov.est    <- tcrossprod(post.load) + diag(post.psi)
        } else {
          cov.est    <- diag(post.psi)
        }
        dimnames(cov.est)      <- list(varnames, varnames)
      } else {
        if(all(!sw["l.sw"], Qg > 0, !sw["psi.sw"]))   {
                                  warning("Loadings & Uniquenesses not stored: can't estimate Sigma and compute error metrics", call.=F)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't estimate Sigma and compute error metrics", call.=F)
        } else if(!sw["psi.sw"])  warning("Uniquenesses not stored: can't estimate Sigma and compute error metrics", call.=F)
      }  
    } else {
      cov.est    <- sims[[G.ind]][[Q.ind]]$cov.est
      if(all(recomp, sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        cov.est  <- replace(cov.est, is.numeric(cov.est), 0)
        for(r in seq_len(n.store))  {
         Sigma   <- tcrossprod(lmat[,,r]) + diag(psi[,r])
         cov.est <- cov.est + Sigma/n.store
        }
      } else if(recomp)   {
        if(all(!sw["l.sw"], Qg > 0, !sw["psi.sw"]))   {
                                  warning("Loadings & Uniquenesses not stored: can't re-estimate Sigma", call.=F)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't re-estimate Sigma", call.=F)
        } else if(!sw["psi.sw"])  warning("Uniquenesses not stored: can't re-estimate Sigma", call.=F)
        
      }
    }
    
    emp.T        <- exists("cov.emp", envir=environment())
    est.T        <- exists("cov.est", envir=environment())
    if(all(emp.T, est.T)) {
      error      <- cov.emp - cov.est
      MSE[g]     <- mean(error * error)
      RMSE[g]    <- sqrt(MSE[g])
      NRMSE[g]   <- RMSE[g]/(max(cov.emp) - min(cov.emp))
      CVRMSE[g]  <- RMSE[g]/mean(cov.emp)
      MAD[g]     <- mean(abs(error))
      if(any(all(isTRUE(scaling), cent)    && 
                 sum(round(diag(cov.est))  != 
                 round(diag(cov.emp)))     != 0,
         sum(abs(post.psi - (1 - post.psi)) < 0) != 0,
         var.exp  > 1))           warning(paste0(ifelse(G == 1, "C", paste0("Group ", g, "'s c")), "hain may not have fully converged"), call.=F)
    }
  
    results      <- list(if(sw["mu.sw"])   list(means = mu), 
                         list(post.mu    = post.mu),
                         if(sw["psi.sw"])  list(uniquenesses = psi), 
                         list(post.psi   = post.psi),
                         if(sw["l.sw"])    list(loadings     = lmat, 
                                                post.load    = post.load),
                         list(var.exp    = var.exp),
                         if(emp.T) list(cov.mat = cov.emp), 
                         if(est.T) list(cov.est = cov.est))
    result[[g]]  <- unlist(results, recursive=F)
    attr(result[[g]], "Store") <- n.store
  }
  names(result)  <- paste0("Group", seq_len(G))
  attr(GQ.res, "Criterion")    <- criterion
  attr(GQ.res, "Factors")      <- n.fac
  attr(GQ.res, "Groups")       <- n.grp
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  if(all(emp.T, est.T)) {
    errors       <- list(MSE = mean(MSE), RMSE = mean(RMSE), NRMSE = mean(NRMSE),
                         CVRMSE = mean(CVRMSE), MAD = mean(MAD))  
  }
  
  result         <- c(result, if(exists("cluster", envir=environment())) list(Clust = cluster), 
                      if(all(emp.T, est.T)) list(Error = errors), list(GQ.results = GQ.res), 
                      if(sw["f.sw"]) list(Scores = scores))
  class(result)                <- "IMIFA"
  attr(result, "Method")       <- attr(sims, "Method")
  attr(result, "Obs")          <- n.obs
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}