#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.imifa       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL, Q = NULL, Q.meth = c("Mode", "Median"), G.meth = c("Mode", "Median"),
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
  n.store        <- length(store)
  tmp.store      <- store
  label.switch   <- attr(sims, "Label.Switch")
  method         <- attr(sims, "Method")
  inf.G          <- is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))
  inf.Q          <- !is.element(method, c("FA", "MFA", "OMFA", "IMFA"))
  n.fac          <- attr(sims, "Factors")
  n.grp          <- attr(sims, "Groups")
  n.obs          <- attr(sims, "Obs")
  n.var          <- attr(sims, "Vars")
  sw             <- attr(sims, "Switch")
  cent           <- attr(sims, "Center")
  scaling        <- attr(sims, "Scaling")
  scal.meth      <- attr(scaling, "Method")
  conf.level     <- as.numeric(conf.level)
  if(conf.level   < 0  && 
     conf.level   > 1)            stop("'conf.level' must be a single number between 0 and 1")
  conf.levels    <- c((1 - conf.level)/2, 1 - (1 - conf.level)/2)
  criterion      <- match.arg(criterion)
  if(all(!is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
     !is.element(criterion, 
     c("aicm", "bicm"))))         stop(paste0("'criterion' should be one of 'aicm' or 'bicm' for the ", method, "method"))
  if(!is.logical(recomp))         stop("'recomp' must be TRUE or FALSE")
  if(any(burnin   > 0, 
     thinning     > 1)) recomp <- T
  
  G.T            <- !missing(G)
  Q.T            <- !missing(Q)
  G.ind          <- 1
  Q.ind          <- 1
  if(inf.G) {
    GQs          <- length(sims[[G.ind]])
    G.store      <- matrix(unlist(lapply(seq_len(GQs), function(gq) sims[[G.ind]][[gq]]$G.store[store])), nr=GQs, nc=length(store))
    G.meth       <- ifelse(missing(G.meth), "Mode", match.arg(G.meth))
  }
  if(G.T)   {
    if(!inf.G) {
      if(!is.element(method, c("FA", "IFA"))) {
        G.ind    <- which(n.grp == G)
      } else if(G > 1)            message(paste0("Forced G=1 for the ", method, " method"))
    } else if(!is.element(G, 
              unique(G.store)))   stop("This G value was not visited during simulation")
    if(all(is.element(method, c("MFA", "MIFA")),
       !is.element(G, n.grp)))    stop("This G value was not used during simulation")
  }
  G              <- ifelse(all(G.T, !is.element(method, c("FA", "IFA"))), G, 1)
  if(Q.T)   {
    if(G.T) {
      if(length(Q) == 1)     Q <- rep(Q, G)
      if(length(Q) != G)          stop(paste0("'Q' must be supplied for each of the ", G, " groups"))
    } else if(length(Q) != 1)     stop("'Q' must be a scalar if G=1 or G is not suppplied")
    if(all(is.element(method, c("FA", "MFA")), !G.T)) {
      Q.ind      <- which(n.fac == Q)
    }
    if(all(is.element(method, c("FA", "MFA")), 
       !is.element(Q, n.fac)))    stop("This Q value was not used during simulation")
    if(all(inf.Q, 
      (Q * (n.fac - Q)) < 0))     stop(paste0("Q can't be greater than the number of factors in ", match.call()$sims))
  } 
  
  if(inf.G) {
    G.tab        <- if(GQs > 1) lapply(apply(G.store, 1, function(x) list(table(x, dnn=NULL))), "[[", 1) else table(G.store, dnn=NULL)
    G.prob       <- if(GQs > 1) lapply(G.tab, prop.table) else prop.table(G.tab)
    G.mode       <- if(GQs > 1) unlist(lapply(G.tab, function(gt) as.numeric(names(gt[gt == max(gt)])[1]))) else as.numeric(names(G.tab[G.tab == max(G.tab)])[1])
    G.med        <- if(GQs > 1) ceiling(apply(G.store, 1, median) * 2)/2 else ceiling(median(G.store) * 2)/2
    if(!G.T) {
      G          <- if(G.meth == "Mode") G.mode else floor(G.med)
    }
    G.CI         <- if(GQs > 1) apply(G.store, 1, function(gs) round(quantile(gs, conf.levels))) else round(quantile(G.store, conf.levels))
    tmp.store    <- if(GQs > 1) lapply(seq_len(GQs), function(gq) store[which(G.store[gq,] == G[gq])]) else store[which(G.store == G)]
    GQ.temp1     <- list(G = G, G.Mode = G.mode, G.Median = G.med, 
                         G.CI = G.CI, G.Probs = G.prob, G.Counts = G.tab)
  }
  
  G.range        <- ifelse(G.T, 1, length(n.grp))
  Q.range        <- ifelse(any(Q.T, all(!is.element(method, c("OMFA", "IMFA")), inf.Q)), 1, length(n.fac))
  crit.mat       <- matrix(NA, nr=G.range, nc=Q.range)
    
  # Retrieve log-likelihoods and/or tune G &/or Q according to criterion
    if(all(G.T, Q.T)) {
      dimnames(crit.mat) <- list(paste0("G", G), paste0("Q", ifelse(length(Q) == 1, Q, "IFA")))
    } else if(G.T)    {
      dimnames(crit.mat) <- list(paste0("G", G), paste0("Q", n.fac))
    } else if(Q.T)    {
      dimnames(crit.mat) <- list(paste0("G", n.grp), paste0("Q", Q))
    } else {
      dimnames(crit.mat) <- list(paste0("G", n.grp), paste0("Q", n.fac))
    }
    if(all(!is.element(method, c("OMFA", "IMFA")),
       inf.Q)) {
      colnames(crit.mat) <- "IFA"
    }
    if(inf.G)  {
      rownames(crit.mat) <- "IM"
    }
    aicm         <- bicm       <- 
    aic.mcmc     <- bic.mcmc   <- crit.mat
    log.N        <- log(n.obs)
    for(g in seq_len(G.range))   { 
      gi                 <- ifelse(G.T, G.ind, g)
      for(q in seq_len(Q.range)) {
        qi               <- ifelse(Q.T, Q.ind, q)
        log.likes        <- if(is.element(method, c("OMFA", "IMFA")) && GQs > 1) sims[[gi]][[qi]]$ll.store[tmp.store[[qi]]] else sims[[gi]][[qi]]$ll.store[tmp.store]
        ll.max           <- 2 * max(log.likes, na.rm=T)
        ll.var           <- ifelse(length(log.likes) != 1, 2 * var(log.likes, na.rm=T), 0)
        ll.mean          <- mean(log.likes, na.rm=T)
        aicm[g,q]        <- ll.max - ll.var * 2
        bicm[g,q]        <- ll.max - ll.var * log.N
        if(!inf.Q) {
          K              <- if(!is.element(method, c("OMFA", "IMFA"))) attr(sims[[gi]][[qi]], "K") else G[qi] - 1 + G[qi] * (n.var * n.fac[qi] - 0.5 * n.fac[qi] * (n.fac[qi] - 1)) + 2 * G[qi] * n.var
          aic.mcmc[g,q]  <- ll.max - K * 2
          bic.mcmc[g,q]  <- ll.max - K * log.N
        }
      }  
    }
    crit         <- get(criterion)
    crit.max     <- which(crit == max(crit), arr.ind = T)
  
  # Control for supplied values of G &/or Q
    if(!any(Q.T, G.T)) {
      G.ind      <- crit.max[1]
      Q.ind      <- crit.max[2]
      if(!inf.G) {
        G        <- n.grp[G.ind]
      }
      if(!inf.Q) {
        Q        <- n.fac[Q.ind]  
      }
    } else if(all(G.T, !Q.T)) {
      Q.ind      <- which(crit == max(crit))
      if(!inf.Q) {
        Q        <- n.fac[Q.ind]  
      }
    } else if(all(Q.T, !G.T)) {
      G.ind      <- which(crit == max(crit))
      if(!inf.G) {
        G        <- n.grp[G.ind]
      }
    } 
    G            <- ifelse(all(length(n.grp) == 1, !inf.G), n.grp, G)
    Gseq         <- seq_len(G)
    Gseq2        <- if(inf.G) seq_len(n.grp) else Gseq
    G.ind        <- ifelse(all(length(n.grp) == 1, !inf.G), which(n.grp == G), G.ind)
    GQ.temp2     <- list(AICMs = aicm, BICMs = bicm)
    if(is.element(method, c("OMFA", "IMFA")) &&
       GQs > 1)  {
      tmp.store  <- tmp.store[[Q.ind]]
    }
    if(!inf.Q)   {
      Q          <- if(length(n.fac) > 1) Q     else n.fac
      Q.ind      <- if(length(n.fac) > 1) Q.ind else which(n.fac == Q)
      Q          <- setNames(rep(Q, G), paste0("Group ", Gseq))  
      GQ.temp1   <- if(is.element(method, c("OMFA", "IMFA")) && GQs > 1) lapply(GQ.temp1, "[[", 1) else if(inf.G) GQ.temp1
      GQ.temp3   <- c(GQ.temp2, list(AIC.mcmcs = aic.mcmc, BIC.mcmcs = bic.mcmc))
      GQ.res     <- if(!is.element(method, c("OMFA", "IMFA"))) c(list(G = G, Q = Q), GQ.temp3) else c(GQ.temp1, list(Q = Q), GQ.temp3)
    }
    if(inf.G)    {
      non.empty  <- sims[[G.ind]][[Q.ind]]$nonempty[tmp.store]
    }
    clust.ind    <- !any(is.element(method,  c("FA", "IFA")), 
                     all(is.element(method, c("MFA", "MIFA")), G == 1))
    sw.mx        <- ifelse(clust.ind, sw["mu.sw"], T)
    sw.px        <- ifelse(clust.ind, sw["psi.sw"], T)  
    if(inf.Q) {
      Q.store    <- sims[[G.ind]][[Q.ind]]$Q.store[,tmp.store, drop=F]
      Q.meth     <- ifelse(missing(Q.meth), "Mode", match.arg(Q.meth))
    }
  
# Manage Label Switching & retrieve cluster labels/mixing proportions
  if(clust.ind) {
    source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=T)
    if(sw["mu.sw"])   {
      mus        <- sims[[G.ind]][[Q.ind]]$mu[,,tmp.store, drop=F]                            
    }
    if(sw["l.sw"])    {
      lmats      <- sims[[G.ind]][[Q.ind]]$load
      if(!is.element(method, c("MFA", "OMFA", "IMFA"))) {
        lmats    <- as.array(lmats)
      }
      lmats      <- lmats[,,,tmp.store, drop=F]
    }
    if(sw["psi.sw"])  {
      psis       <- sims[[G.ind]][[Q.ind]]$psi[,,tmp.store, drop=F]
    }
    if(sw["pi.sw"])   {
      pies       <- sims[[G.ind]][[Q.ind]]$pi.prop[,tmp.store, drop=F]
    }
    z            <- as.matrix(sims[[G.ind]][[Q.ind]]$z.store[,tmp.store])
    label.miss   <- missing(Labels)
    if(!label.miss)   {
      if(!exists(as.character(substitute(Labels)),
          envir=.GlobalEnv))      stop(paste0("Object ", match.call()$Labels, " not found"))
      if(length(Labels) != n.obs) stop(paste0("'Labels' must be a factor of length N=",  n.obs))
      zlabels    <- factor(Labels, labels=seq_along(unique(Labels)))
      levs       <- levels(zlabels)
      lev.ind    <- length(levs) == G
    }
    if(!label.switch) {
      z.temp     <- factor(z[,1], labels=Gseq)
      if(!label.miss && lev.ind) {    
        sw.lab   <- lab.switch(z.new=z.temp, z.old=zlabels, Gs=Gseq2)
        z.temp   <- sw.lab$z
        l.perm   <- sw.lab$z.perm
      }
      for(sl in seq_along(tmp.store)[-1])   {
        n.ind    <- if(inf.G) non.empty[[sl]] else Gseq2
        Nseq     <- seq_along(n.ind)
        sw.lab   <- lab.switch(z.new=z[,sl], z.old=z.temp, Gs=Nseq)
        z[,sl]   <- sw.lab$z
        z.perm   <- sw.lab$z.perm
        perm     <- identical(as.integer(z.perm), Nseq)
        if(!perm) {
          if(sw["mu.sw"])  {
            mus[,,sl]    <- mus[,z.perm,sl]
          }
          if(sw["l.sw"])   {
            lmats[,,,sl] <- lmats[,,z.perm,sl]
          }
          if(sw["psi.sw"]) {
            psis[,,sl]   <- psis[,z.perm,sl]
          }
          if(sw["pi.sw"])  {
            pies[,sl]    <- pies[z.perm,sl]
          }
          if(inf.Q)        {
            Q.store[,sl] <- Q.store[z.perm,sl]
          }  
        }
      }
    }
    post.z       <- droplevels(setNames(apply(z, 1, function(x) factor(which.max(tabulate(x)), levels=Gseq2)), seq_len(n.obs)))
    var.z        <- apply(z, 1, var)
    CI.z         <- apply(z, 1, function(x) round(quantile(x, conf.levels)))
    if(sw["pi.sw"])    {
      pi.prop    <- pies[,seq_along(tmp.store)]
      post.pi    <- rowMeans(pi.prop, dims=1)
      var.pi     <- apply(pi.prop, 1, var)
      CI.pi      <- apply(pi.prop, 1, function(x) quantile(x, conf.levels))
    } else {
      post.pi    <- setNames(prop.table(tabulate(post.z, nbins=G)), paste0("Group ", Gseq))
    }
    if(!label.miss) {
      if(nlevels(post.z) == length(levs)) {
        sw.lab   <- lab.switch(z.new=post.z, z.old=zlabels, Gs=Gseq2)
        post.z   <- factor(sw.lab$z)
        l.perm   <- sw.lab$z.perm
        post.pi  <- post.pi[l.perm]  
      }
      tab        <- table(post.z, zlabels, dnn=list("Predicted", "Observed"))
      tab.stat   <- classAgreement(tab)
    }
    cluster      <- list(post.z = post.z, post.pi = post.pi, 
                         z = z, var.z = var.z, CI.z = CI.z)
    cluster      <- c(cluster, if(!label.miss) list(conf.mat = tab, perf = tab.stat),
                      if(sw["pi.sw"]) list(pi.prop = pi.prop, var.pi = var.pi, CI.pi = CI.pi))
    attr(cluster, "Z.init")    <- attr(sims[[G.ind]], "Z.init")
    attr(cluster, "Init.Meth") <- attr(sims, "Init.Z")
    post.z       <- as.numeric(levels(post.z))[post.z]
    z.ind        <- lapply(Gseq, function(g) post.z == g)
  }
  if(inf.Q)   {
    Q.store      <- if(G > 1) Q.store[Gseq,] else Q.store
    Q.tab        <- if(G > 1) lapply(apply(Q.store, 1, function(x) list(table(x, dnn=NULL))), "[[", 1) else table(Q.store, dnn=NULL)
    Q.prob       <- if(G > 1) lapply(Q.tab, prop.table) else prop.table(Q.tab)
    Q.mode       <- if(G > 1) unlist(lapply(Q.tab, function(qt) as.numeric(names(qt[qt == max(qt)])[1]))) else as.numeric(names(Q.tab[Q.tab == max(Q.tab)])[1])
    Q.med        <- if(G > 1) ceiling(apply(Q.store, 1, median) * 2)/2 else ceiling(median(Q.store) * 2)/2
    if(!Q.T)  {
      Q          <- if(Q.meth == "Mode") Q.mode else floor(Q.med)
    } else    {
      Q          <- if(G.T) Q else setNames(rep(Q, G), paste0("Group ", Gseq))
    }
    Q.CI         <- if(G > 1) apply(Q.store, 1, function(qs) round(quantile(qs, conf.levels))) else round(quantile(Q.store, conf.levels))
    GQ.temp4     <- list(Q = Q, Q.Mode = Q.mode, Q.Median = Q.med, 
                         Q.CI = Q.CI, Q.Probs = Q.prob, Q.Counts = Q.tab)
    if(inf.G) {
      GQ.res     <- c(GQ.temp1, GQ.temp4)
    } else    {
      GQ.res     <- c(list(G = G), GQ.temp4)
    }
    GQ.res       <- c(GQ.res, GQ.temp2)
  }

# Retrieve (unrotated) scores
  if(all(Q == 0)) {
    if(sw["f.sw"])                warning("Scores & loadings not stored as model has zero factors", call.=F)
    sw["f.sw"]   <- F
    no.score     <- T
  }
  if(sw["f.sw"])  {
    Q.max        <- max(Q) 
    Q.maxs       <- seq_len(Q.max)
    f            <- sims[[G.ind]][[Q.ind]]$f
    if(inf.Q) {
      f          <- as.array(f)
    }
    f            <- f[,Q.maxs,tmp.store, drop=F]
  }

# Loop over g in G to extract other results
  result         <- list(list())
  f.store        <- list()
  MSE   <- RMSE  <- NRMSE  <- CVRMSE  <- 
  MAD   <- emp.T <- est.T  <- rep(NA, G)
  for(g in Gseq) {
    Qg           <- Q[g]
    Qgs          <- seq_len(Qg)
    sw["l.sw"]   <- attr(sims, "Switch")["l.sw"]
    if(Qg == 0)  {
      if(all(sw["l.sw"],
             !no.score))          warning(paste0("Loadings not stored as", ifelse(G > 1, paste0(" group ", g), " model"), " has zero factors"), call.=F)
      sw["l.sw"] <- F
    }
    if(inf.Q) {
      store      <- seq_along(tmp.store)[which(Q.store[g,] >= Qg)]
    } else    {
      store      <- seq_along(tmp.store)
    }
    n.store      <- length(store)
  
  # Retrieve (unrotated) loadings  
    if(sw["l.sw"]) {
      if(clust.ind)  {
        lmat     <- adrop(lmats[,Qgs,g,store, drop=F], drop=3)
        l.temp   <- adrop(lmat[,,1, drop=F], drop=3)
      }
      if(!clust.ind) {
        lmat     <- sims[[G.ind]][[Q.ind]]$load
        if(inf.Q) {
          lmat   <- as.array(lmat)
        }
        lmat     <- lmat[,Qgs,store, drop=F]
        l.temp   <- adrop(lmat[,,1, drop=F], drop=3)
      }
    }
  
  # Loadings matrix / identifiability / error metrics / etc.  
    if(all(sw["f.sw"], clust.ind)) {
      fg         <- f[z.ind[[g]],Qgs,, drop=F]
    }
    if(sw["l.sw"])      {
      for(p in seq_len(n.store)) {
        pf               <- store[p]
        rot              <- procrustes(X=as.matrix(lmat[,,p]), Xstar=l.temp)$R
        lmat[,,p]        <- lmat[,,p] %*% rot
        if(sw["f.sw"])  {
          if(clust.ind) {
            fg[,,pf]     <- fg[,,pf]  %*% rot
          } else {
            f[,,pf]      <- f[,,pf]   %*% rot
          }
        }  
      }
    }
    if(all(sw["f.sw"], clust.ind)) {
      f[z.ind[[g]],Qgs,] <- fg
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
        dat      <- attr(sims, "Name")
        data.x   <- exists(dat, envir=.GlobalEnv)  
        if(!data.x)               warning(paste0("Object ", dat, " not found in .GlobalEnv: can't compute empirical covariance and error metrics"), call.=F) 
        dat      <- as.data.frame(get(dat))
        dat      <- dat[sapply(dat, is.numeric)]
        dat      <- scale(dat, center=cent, scale=scaling)
        varnames <- colnames(dat)
      }
      cov.emp    <- cov(dat[z.ind[[g]],, drop=F])
      dimnames(cov.emp)  <- list(varnames, varnames)
      if(sum(z.ind[[g]]) == 0)    rm(cov.emp)
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
      class(post.load)   <- "loadings"
    } else   {
      var.exp    <- ifelse(sum(z.ind[[g]]) == 0, 0, max(0, (sum(diag(cov.emp)) - sum(post.psi))/n.var))
    }
  
  # Calculate estimated covariance matrices & compute error metrics
    if(clust.ind) {
      if(all(sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        if(Qg > 0)   {
          cov.est        <- tcrossprod(post.load) + diag(post.psi)
        } else {
          cov.est        <- diag(post.psi)
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
    f.store[[g]] <- store
  }
  if(sw["f.sw"]) {
    f            <- f[,,unique(unlist(f.store)), drop=F]
    scores       <- list(f = f, post.f = rowMeans(f, dims=2), var.f = apply(f, c(1, 2), var),
                         CI.f  = apply(f, c(1, 2), function(x) quantile(x, conf.levels)))
  }
  names(result)  <- paste0("Group", Gseq)
  class(GQ.res)                <- "listof"
  attr(GQ.res, "Criterion")    <- criterion
  attr(GQ.res, "Factors")      <- n.fac
  attr(GQ.res, "Groups")       <- n.grp
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  err.T                        <- unlist(lapply(Gseq, function(g) all(emp.T[g], est.T[g])))
  if(any(err.T)) {
    errors       <- list(MSE = mean(MSE, na.rm=T), RMSE = mean(RMSE, na.rm=T), NRMSE = mean(NRMSE, na.rm=T),
                         CVRMSE = mean(CVRMSE, na.rm=T), MAD = mean(MAD, na.rm=T))  
    class(errors)              <- "listof"
  }
  
  result         <- c(result, if(exists("cluster", envir=environment())) list(Clust = cluster), 
                      if(any(err.T)) list(Error = errors), list(GQ.results = GQ.res), 
                      if(sw["f.sw"]) list(Scores = scores))
  class(result)                <- "IMIFA"
  attr(result, "Method")       <- method
  attr(result, "Name")         <- attr(sims, "Name")
  attr(result, "Obs")          <- n.obs
  attr(result, "Store")        <- tmp.store
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}