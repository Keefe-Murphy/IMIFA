#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.imifa       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL, Q = NULL, Q.meth = c("Mode", "Median"),
                             criterion = c("bicm", "aicm", "bic.mcmc", "aic.mcmc"), conf.level = 0.95, Labels = NULL, recomp = F) {
  
  defpar         <- suppressWarnings(par(no.readonly = T))
  defop          <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(missing(sims))               stop("Simulations must be supplied")
  if(!exists(deparse(substitute(sims)),
             envir=.GlobalEnv))   stop(paste0("Object ", match.call()$sims, " not found"))
  if(class(sims) != "IMIFA")      stop(paste0("Simulations object of class 'IMIFA' must be supplied"))
  burnin         <- as.integer(burnin)
  thinning       <- as.integer(thinning)
  store          <- seq(from=burnin + 1, to=attr(sims, "Store"), by=thinning)
  temp.store     <- store
  n.store        <- length(store)
  method         <- attr(sims, "Method")
  n.fac          <- attr(sims, "Factors")
  n.grp          <- attr(sims, "Groups")
  n.obs          <- attr(sims, "Obs")
  n.var          <- attr(sims, "Vars")
  sw             <- attr(sims, "Switch")
  cent           <- attr(sims, "Center")
  scaling        <- attr(sims, "Scaling")
  scal.meth      <- attr(scaling, "Method")
  conf.level     <- as.numeric(conf.level)
  if(abs(conf.level -
        (1 - conf.level)) < 0)    stop("'conf.level' must be a single number between 0 and 1")
  conf.levels    <- c((1 - conf.level)/2, 1 - (1 - conf.level)/2)
  criterion      <- match.arg(criterion)
  if(all(method  == "MIFA", 
     !is.element(criterion, 
     c("aicm", "bicm"))))         stop("'criterion' should be one of 'aicm' or 'bicm' for the MIFA method")
  if(!is.logical(recomp))         stop("'recomp' must be TRUE or FALSE")
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
      (Q * (n.fac - Q)) < 0))     stop(paste0("Q can't be greater than the number of factors in ", match.call()$sims))
  } 
  G              <- ifelse(all(G.T, !is.element(method, c("FA", "IFA"))), G, 1)
  
  if(is.element(method, c("IFA", "MIFA"))) {
    if(missing(Q.meth)) {
      Q.meth     <- "Mode"
    } else   {
      Q.meth     <- match.arg(Q.meth)
    }
    
  # Set Q as the (lesser of) the distribution's mode(s) & compute credible interval
    Q.store      <- sims[[G.ind]][[Q.ind]]$Q.store[,store, drop=F]
    NQ           <- nrow(Q.store)
    Q.tab        <- if(NQ > 1) lapply(apply(Q.store, 1, function(x) list(table(x, dnn=NULL))), "[[", 1) else table(Q.store, dnn=NULL)
    Q.prob       <- if(NQ > 1) lapply(Q.tab, prop.table) else prop.table(Q.tab)
    Q.mode       <- if(NQ > 1) unlist(lapply(Q.tab, function(qt) as.numeric(names(qt[qt == max(qt)])[1]))) else as.numeric(names(Q.tab[Q.tab == max(Q.tab)])[1])
    Q.med        <- ceiling(apply(Q.store, 1, median) * 2)/2
    Q            <- if(Q.meth == "Mode") Q.mode else Q.med
    Q.CI         <- if(NQ > 1) apply(Q.store, 1, function(qs) round(quantile(qs, conf.levels))) else round(quantile(Q.store, conf.levels))
    GQ.res       <- list(G = G, Q = Q, Mode = Q.mode, Median = Q.med, 
                         Q.CI = Q.CI, Probs= Q.prob, Counts = Q.tab)
    if(method == "MIFA") {
      GQres.temp <- GQ.res[-seq_len(2)]
    }
  }
  
  if(is.element(method, c("FA", "MFA", "MIFA"))) {
    G.range      <- ifelse(G.T, 1, length(n.grp))
    Q.range      <- ifelse(any(Q.T, method == "MIFA"), 1, length(n.fac))
    crit.mat     <- matrix(NA, nr=G.range, nc=Q.range)
    
  # Retrieve log-likelihoods and tune G & Q according to criterion
    if(all(G.T, Q.T)) {
      dimnames(crit.mat) <- list(paste0("G", G), paste0("Q", Q))
    } else if(G.T)    {
      dimnames(crit.mat) <- list(paste0("G", G), paste0("Q", n.fac))
    } else if(Q.T)    {
      dimnames(crit.mat) <- list(paste0("G", n.grp), paste0("Q", Q))
    } else {
      dimnames(crit.mat) <- list(paste0("G", n.grp), paste0("Q", n.fac))
    }
    if(method == "MIFA") {
      colnames(crit.mat) <- "IFA"
    }
    aicm         <- bicm       <- 
    aic.mcmc     <- bic.mcmc   <- crit.mat
    for(g in seq_len(G.range)) { 
      gi                 <- ifelse(G.T, G.ind, g)
      for(q in seq_len(Q.range)) {
        qi               <- ifelse(Q.T, Q.ind, q)
        log.likes        <- sims[[gi]][[qi]]$ll.store[store]
        ll.max           <- 2 * max(log.likes, na.rm=T)
        ll.var           <- ifelse(length(log.likes) != 1, 2 * var(log.likes, na.rm=T), 0)
        ll.mean          <- mean(log.likes, na.rm=T)
        aicm[g,q]        <- ll.max - ll.var * 2
        bicm[g,q]        <- ll.max - ll.var * log(n.obs) 
        if(method != "MIFA") {
          K              <- attr(sims[[gi]][[qi]], "K")
          aic.mcmc[g,q]  <- ll.max - K * 2
          bic.mcmc[g,q]  <- ll.max - K * log(n.obs)
        }
      }  
    }
    crit         <- get(criterion)
    crit.max     <- which(crit == max(crit), arr.ind = T)
  
  # Control for supplied values of G &/or Q
    if(!any(Q.T, G.T)) {
      G.ind      <- crit.max[1]
      Q.ind      <- crit.max[2]
      G          <- n.grp[G.ind]
      Q          <- if(method == "MIFA") Q else n.fac[Q.ind]
    } else if(all(G.T, !Q.T)) {
      Q.ind      <- which(crit == max(crit))
      Q          <- if(method == "MIFA") Q else n.fac[Q.ind]
    } else if(all(Q.T, !G.T)) {
      G.ind      <- which(crit == max(crit))
      G          <- n.grp[G.ind]
    } 
    G            <- ifelse(length(n.grp) == 1, n.grp, G)
    Q            <- if(any(method == "MIFA", length(n.fac) > 1)) Q else n.fac
    G.ind        <- ifelse(length(n.grp) == 1, which(n.grp == G), G.ind)
    Q.ind        <- if(any(method == "MIFA", length(n.fac) > 1)) Q.ind else which(n.fac == Q)
    Q            <- if(method == "MIFA") Q else setNames(rep(Q, G), paste0("Group ", seq_len(G)))
    GQ.res       <- list(G = G, Q = Q, AICMs = aicm, BICMs = bicm)
    if(method == "MIFA") {
      GQ.res     <- c(GQ.res, GQres.temp)
    } else {
      GQ.res     <- c(GQ.res, list(AIC.mcmcs = aic.mcmc, BIC.mcmcs = bic.mcmc))
    }
  }
  clust.ind      <- all(is.element(method, c("MFA", "MIFA", "IMIFA")), G > 1)
  sw.mx          <- ifelse(clust.ind, sw["mu.sw"], T)
  sw.px          <- ifelse(clust.ind, sw["psi.sw"], T)
  
# Manage Label Switching & retrieve cluster labels/mixing proportions
  if(clust.ind) {
    if(sw["mu.sw"])  {
      mus        <- sims[[G.ind]][[Q.ind]]$mu[,,store, drop=F]                            
    }
    if(sw["l.sw"])   {
      lmats      <- sims[[G.ind]][[Q.ind]]$load
      if(method  == "MIFA") {
        lmats    <- as.array(lmats)
      }
      lmats      <- lmats[,,,store, drop=F]
    }
    if(sw["psi.sw"]) {
      psis       <- sims[[G.ind]][[Q.ind]]$psi[,,store, drop=F]
    }
    if(sw["pi.sw"])  {
      pies       <- sims[[G.ind]][[Q.ind]]$pi.prop[,store, drop=F]
    }
    z            <- as.matrix(sims[[G.ind]][[Q.ind]]$z.store[,store])
    z.temp       <- factor(z[,1], levels=seq_len(G))
    for(ls in seq_len(n.store)[-1]) {
      tab        <- table(factor(z[,ls], levels=seq_len(G)), z.temp)
      z.perm     <- matchClasses(tab, method="exact", verbose=F)
      z[,ls]     <- factor(z[,ls], labels=z.perm, levels=seq_len(G))
      if(sw["mu.sw"])  {
        mus[,,ls]      <- mus[,z.perm,ls]
      }
      if(sw["l.sw"])   {
        lmats[,,,ls]   <- lmats[,,z.perm,ls]
      }
      if(sw["psi.sw"]) {
        psis[,,ls]     <- psis[,z.perm,ls]
      }
      if(sw["pi.sw"])  {
        pies[,ls]      <- pies[z.perm,ls]
      }
    }
    post.z       <- setNames(apply(z, 1, function(x) factor(which.max(tabulate(x)), levels=seq_len(G))), seq_len(n.obs))
    var.z        <- apply(z, 1, var)
    CI.z         <- apply(z, 1, function(x) round(quantile(x, conf.levels)))
    if(sw["pi.sw"])    {
      pi.prop    <- pies[,store]
      post.pi    <- rowMeans(pi.prop, dims=1)
      var.pi     <- apply(pi.prop, 1, var)
      CI.pi      <- apply(pi.prop, 1, function(x) quantile(x, conf.levels))
    } else {
      post.pi    <- setNames(prop.table(tabulate(post.z, nbins=G)), paste0("Group ", seq_len(G)))
    }
    if(!missing(Labels)) {
      if(!exists(as.character(substitute(Labels)),
         envir=.GlobalEnv))           stop(paste0("Object ", match.call()$Labels, " not found"))
      labels     <- as.factor(Labels)
      levs       <- levels(labels)
      if(length(labels)  != n.obs)    stop(paste0("Labels must be a factor of length N=",  n.obs))
      tab        <- table(post.z, as.numeric(labels))
      if(nlevels(post.z) == length(levs)) {
        perm     <- matchClasses(tab, method="exact", verbose=F)
        post.nn  <- tabulate(post.z, nbins=G)
        post.z   <- factor(factor(post.z, labels=levs[perm][post.nn > 0]), levels=levs)
        post.pi  <- post.pi[perm]  
      }
      tab        <- table(post.z, labels, dnn=list("Predicted", "Observed"))
      tab.stat   <- classAgreement(tab)
    }
    cluster      <- list(post.z = post.z, post.pi = post.pi, 
                         z = z, var.z = var.z, CI.z = CI.z)
    cluster      <- c(cluster, if(!missing(Labels)) list(conf.mat = tab, perf = tab.stat),
                      if(sw["pi.sw"]) list(pi.prop = pi.prop, var.pi = var.pi, CI.pi = CI.pi))
    attr(cluster, "Z.init")    <- attr(sim[[G.ind]], "Z.init")
    attr(cluster, "Init.Meth") <- attr(sims, "Init.Z")
    post.z       <- as.numeric(post.z)
    ind          <- lapply(seq_len(G), function(g) post.z == g)
  }
  
# Retrieve (unrotated) scores
  if(all(Q == 0)) {
    if(sw["f.sw"])                warning("Scores not stored as model has zero factors", call.=F)
    sw["f.sw"]   <- F
  }
  if(sw["f.sw"])  {
    Q.max        <- max(Q) 
    Q.maxs       <- seq_len(Q.max)
    f            <- sims[[G.ind]][[Q.ind]]$f
    if(is.element(method, c("IFA", "MIFA"))) {
      f          <- as.array(f)
    }
    f            <- f[,Q.maxs,store, drop=F]
  }

# Loop over g in G to extract other results
  result         <- list(list())
  MSE   <- RMSE  <- NRMSE  <- CVRMSE  <- 
  MAD   <- emp.T <- est.T  <- rep(NA, G)
  for(g in seq_len(G)) {
    Qg           <- Q[g]
    Qgs          <- seq_len(Qg)
    sw["l.sw"]   <- attr(sims, "Switch")["l.sw"]
    if(Qg == 0) {
      if(sw["l.sw"])              warning(paste0("Loadings not stored as", ifelse(G > 1, paste0(" group ", g), " model"), " has zero factors"), call.=F)
      sw["l.sw"] <- F
    }
    if(is.element(method, c("IFA", "MIFA"))) {
      store      <- temp.store[which(Q.store[g,] >= Qg)]
      n.store    <- length(store)
    }
  
  # Retrieve (unrotated) loadings  
    if(sw["l.sw"]) {
      if(all(is.element(method, c("MFA", "MIFA")), G > 1)) {
        lmat     <- adrop(lmats[,Qgs,g,store, drop=F], drop=3)
        l.temp   <- adrop(lmat[,,1, drop=F], drop=3)
      }
      if(any(is.element(method, c("FA", "IFA")), 
         all(is.element(method, c("MFA", "MIFA")), G == 1))) {
        lmat     <- sims[[G.ind]][[Q.ind]]$load
        if(is.element(method, c("IFA", "MIFA"))) {
          lmat   <- as.array(lmat)
        }
        lmat     <- lmat[,Qgs,store, drop=F]
        l.temp   <- adrop(lmat[,,1, drop=F], drop=3)
      }
    }
  
  # Loadings matrix / identifiability / error metrics / etc.  
    if(all(sw["f.sw"], clust.ind)) {
      fg         <- f[ind[[g]],Qgs,, drop=F]
    }
    if(sw["l.sw"])     {
      for(p in seq_len(n.store)) {
        rot            <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
        lmat[,,p]      <- lmat[,,p] %*% rot
        if(sw["f.sw"])  {
          if(clust.ind) {
            fg[,,p]    <- fg[,,p]   %*% rot
          } else {
            f[,,p]     <- f[,,p]    %*% rot
          }
        }  
      }
    }
    if(all(sw["f.sw"], clust.ind)) {
      f[ind[[g]],Qgs,] <- fg
    }
  
  # Retrieve means, uniquenesses & empirical covariance matrix
    if(clust.ind) {
      if(sw["mu.sw"])  {
        mu       <- as.matrix(mus[,g,store])
      }
      if(sw["psi.sw"]) {
        psi      <- as.matrix(psis[,g,store])
      }
      if(g == 1) {
        data     <- attr(sims, "Name")
        data.x   <- exists(data, envir=.GlobalEnv)  
        if(!data.x)               warning(paste0("Object ", data, " not found in .GlobalEnv: can't compute empirical covariance and error metrics"), call.=F) 
        data     <- as.data.frame(get(data))
        data     <- data[sapply(data, is.numeric)]
        data     <- scale(data, center=cent, scale=scaling)
        varnames <- colnames(data)
      }
      cov.emp    <- cov(data[ind[[g]],, drop=F])
      dimnames(cov.emp)        <- list(varnames, varnames)
      if(sum(ind[[g]]) == 0)      rm(cov.emp)
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
    if(sw["mu.sw"])  {
      post.mu    <- rowMeans(mu, dims=1)
      var.mu     <- apply(mu, 1, var)
      CI.mu      <- apply(mu, 1, function(x) quantile(x, conf.levels))
    }
    if(sw["psi.sw"]) {
      post.psi   <- rowMeans(psi, dims=1)
      var.psi    <- apply(psi, 1, var)
      CI.psi     <- apply(psi, 1, function(x) quantile(x, conf.levels))
    }
    if(sw["l.sw"])   { 
      post.load  <- rowMeans(lmat, dims=2)
      var.load   <- apply(lmat, c(1, 2), var)
      CI.load    <- apply(lmat, c(1, 2), function(x) quantile(x, conf.levels))
      var.exp    <- sum(colSums(post.load * post.load))/n.var
      class(post.load) <- "loadings"
    } else   {
      var.exp    <- ifelse(sum(z.ind[[g]]) == 0, 0, max(0, (sum(diag(cov.emp)) - sum(post.psi))/n.var))
    }
  
  # Calculate estimated covariance matrices & compute error metrics
    if(clust.ind) {
      if(all(sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        if(Qg > 0)   {
          cov.est      <- tcrossprod(post.load) + diag(post.psi)
        } else {
          cov.est      <- diag(post.psi)
        }
        if(data.x)      {
          dimnames(cov.est)    <- list(varnames, varnames)
        }
      } else if(g == 1) {
        if(all(!sw["l.sw"], Qg  > 0, !sw["psi.sw"]))  {
                                  warning("Loadings & Psi not stored: can't estimate Sigma and compute error metrics", call.=F)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't estimate Sigma and compute error metrics", call.=F)
        } else if(!sw["psi.sw"])  warning("Psi not stored: can't estimate Sigma and compute error metrics", call.=F)
      }  
    } else     {
      cov.est    <- sims[[G.ind]][[Q.ind]]$cov.est
      if(all(recomp, sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        cov.est  <- replace(cov.est, is.numeric(cov.est), 0)
        for(r in seq_len(n.store))    {
          if(Qg > 0) {
            Sig  <- tcrossprod(lmat[,,r]) + diag(psi[,r])
          } else {
            Sig  <- diag(psi[,r])
          }
         cov.est <- cov.est + Sig/n.store
        }
      } else if(all(recomp,  g == 1)) {
        if(all(!sw["l.sw"], Qg  > 0, !sw["psi.sw"]))  {
                                  warning("Loadings & Psi not stored: can't re-estimate Sigma", call.=F)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't re-estimate Sigma", call.=F)
        } else if(!sw["psi.sw"])  warning("Psi not stored: can't re-estimate Sigma", call.=F) 
      }
    }
    
    emp.T[g]     <- exists("cov.emp", envir=environment())
    est.T[g]     <- exists("cov.est", envir=environment())
    if(all(emp.T[g], est.T[g])) {
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
    
    results      <- list(if(sw["mu.sw"])  list(means     = mu,
                                               var.mu    = var.mu,
                                               CI.mu     = CI.mu), 
                         if(sw["l.sw"])   list(loadings  = lmat, 
                                               post.load = post.load,
                                               var.load  = var.load,
                                               CI.load   = CI.load),
                         if(sw["psi.sw"]) list(psi       = psi,
                                               var.psi   = var.psi,
                                               CI.psi    = CI.psi),
                         if(sw.mx)        list(post.mu   = post.mu), 
                         if(sw.px)        list(post.psi  = post.psi),
                         if(any(sw["l.sw"], 
                                sw.px))   list(var.exp   = var.exp),
                         if(emp.T[g])     list(cov.mat   = cov.emp), 
                         if(est.T[g])     list(cov.est   = cov.est))
    result[[g]]  <- unlist(results, recursive=F)
    attr(result[[g]], "Store") <- n.store
  }
  if(sw["f.sw"]) {
    scores       <- list(f = f, post.f = rowMeans(f, dims=2), var.f = apply(f, c(1, 2), var),
                         CI.f  = apply(f, c(1, 2), function(x) quantile(x, conf.levels)))
  }
  names(result)  <- paste0("Group", seq_len(G))
  attr(GQ.res, "Criterion")    <- criterion
  attr(GQ.res, "Factors")      <- n.fac
  attr(GQ.res, "Groups")       <- n.grp
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  err.T                        <- unlist(lapply(seq_len(G), function(g) all(emp.T[g], est.T[g])))
  if(any(err.T)) {
    errors       <- list(MSE = mean(MSE, na.rm=T), RMSE = mean(RMSE, na.rm=T), NRMSE = mean(NRMSE, na.rm=T),
                         CVRMSE = mean(CVRMSE, na.rm=T), MAD = mean(MAD, na.rm=T))  
  }
  
  result         <- c(result, if(exists("cluster", envir=environment())) list(Clust = cluster), 
                      if(any(err.T)) list(Error = errors), list(GQ.results = GQ.res), 
                      if(sw["f.sw"]) list(Scores = scores))
  class(result)                <- "IMIFA"
  attr(result, "Method")       <- attr(sims, "Method")
  attr(result, "Name")         <- attr(sims, "Name")
  attr(result, "Obs")          <- n.obs
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}