#################################################
### Tune Parameters (Single & Shrinkage Case) ###
#################################################

tune.IMIFA       <- function(sims = NULL, burnin = 0, thinning = 1, G = NULL, Q = NULL, Q.meth = c("Mode", "Median"), G.meth = c("Mode", "Median"),
                             criterion = c("bicm", "aicm", "bic.mcmc", "aic.mcmc"), conf.level = 0.95, labels = NULL, recomp = FALSE) {
  
  source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=TRUE)
  defpar         <- suppressWarnings(par(no.readonly=TRUE))
  defopt         <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
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
  MH.step        <- attr(sims, "MH.step")
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
  if(conf.level  <= 0   || 
     conf.level  >= 1)            stop("'conf.level' must be a single number between 0 and 1")
  conf.levels    <- c((1 - conf.level)/2, (1 + conf.level)/2)
  criterion      <- match.arg(criterion)
  if(all(!is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
     !is.element(criterion, 
     c("aicm", "bicm"))))         stop(paste0("'criterion' should be one of 'aicm' or 'bicm' for the ", method, "method"))
  if(!is.logical(recomp))         stop("'recomp' must be TRUE or FALSE")
  if(any(burnin   > 0, 
     thinning     > 1)) recomp <- TRUE
  
  G.T            <- !missing(G)
  Q.T            <- !missing(Q)
  G.ind          <- 1
  Q.ind          <- 1
  if(inf.G) {
    GQs          <- length(sims[[G.ind]])
    G.store      <- matrix(unlist(lapply(seq_len(GQs), function(gq) sims[[G.ind]][[gq]]$G.store[store])), nr=GQs, nc=length(store), byrow=TRUE)
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
      (Q * (n.fac - Q))  < 0))    stop(paste0("Q can't be greater than the number of factors in ", match.call()$sims))
    if(all(inf.Q, 
      (Q * (n.var - Q)) <= 0))    stop(paste0("Q must be less than the number of variables ", n.var))
  } 
  
  if(inf.G)  {
    GQ1          <- GQs  > 1
    G.tab        <- if(GQ1) lapply(apply(G.store, 1, function(x) list(table(x, dnn=NULL))), "[[", 1) else table(G.store, dnn=NULL)
    G.prob       <- if(GQ1) lapply(G.tab, prop.table) else prop.table(G.tab)
    G.mode       <- if(GQ1) unlist(lapply(G.tab, function(gt) as.numeric(names(gt[gt == max(gt)])[1]))) else as.numeric(names(G.tab[G.tab == max(G.tab)])[1])
    G.med        <- if(GQ1) ceiling(apply(G.store, 1, median) * 2)/2 else ceiling(median(G.store) * 2)/2
    if(!G.T) {
      G          <- if(G.meth == "Mode") G.mode else floor(G.med)
    }
    G.CI         <- if(GQ1) apply(G.store, 1, function(gs) round(quantile(gs, conf.levels))) else round(quantile(G.store, conf.levels))
    tmp.store    <- if(GQ1) lapply(seq_len(GQs), function(gq) store[which(G.store[gq,] == G[gq])]) else store[which(G.store == G)]
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
    if(inf.Q) {
      colnames(crit.mat) <- "IFA"
    }
    if(inf.G) {
      rownames(crit.mat) <- "IM"
    }
    aicm         <- bicm       <- 
    aic.mcmc     <- bic.mcmc   <- crit.mat
    log.N        <- log(n.obs)
    for(g in seq_len(G.range))   { 
      gi                 <- ifelse(G.T, G.ind, g)
      for(q in seq_len(Q.range)) {
        qi               <- ifelse(Q.T, Q.ind, q)
        log.likes        <- if(is.element(method, c("OMFA", "IMFA")) && GQ1) sims[[gi]][[qi]]$ll.store[tmp.store[[qi]]] else sims[[gi]][[qi]]$ll.store[tmp.store]
        ll.max           <- 2 * max(log.likes, na.rm=TRUE)
        ll.var           <- ifelse(length(log.likes) != 1, 2 * var(log.likes, na.rm=TRUE), 0)
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
    crit.max     <- which(crit == max(crit), arr.ind=TRUE)
  
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
    G            <- ifelse(inf.G, G[Q.ind], ifelse(length(n.grp) == 1, n.grp, G))
    Gseq         <- seq_len(G)
    G.ind        <- ifelse(all(length(n.grp) == 1, !inf.G), which(n.grp == G), G.ind)
    GQ.temp2     <- list(AICMs = aicm, BICMs = bicm)
    if(is.element(method, c("OMFA", "IMFA")) &&
       GQ1)      {
      tmp.store  <- tmp.store[[Q.ind]]
    }
    if(!inf.Q)   {
      Q          <- if(length(n.fac) > 1) Q     else n.fac
      Q.ind      <- if(length(n.fac) > 1) Q.ind else which(n.fac == Q)
      Q          <- setNames(rep(Q, G), paste0("Group ", Gseq))  
      GQ.temp1   <- if(is.element(method, c("OMFA", "IMFA")) && GQs > 1) lapply(GQ.temp1, "[[", Q.ind) else if(inf.G) GQ.temp1
      GQ.temp3   <- c(GQ.temp2, list(AIC.mcmcs = aic.mcmc, BIC.mcmcs = bic.mcmc))
      GQ.res     <- if(!is.element(method, c("OMFA", "IMFA"))) c(list(G = G, Q = Q), GQ.temp3) else c(GQ.temp1, list(Q = Q), GQ.temp3)
    }
    clust.ind    <- !any(is.element(method,  c("FA", "IFA")), 
                     all(is.element(method, c("MFA", "MIFA")), G == 1))
    sw.mx        <- ifelse(clust.ind, sw["mu.sw"], T)
    sw.px        <- ifelse(clust.ind, sw["psi.sw"], T)  
    if(inf.Q) {
      Q.store    <- sims[[G.ind]][[Q.ind]]$Q.store[Gseq,tmp.store, drop=FALSE]
      Q.meth     <- ifelse(missing(Q.meth), "Mode", match.arg(Q.meth))
    }
  
# Manage Label Switching & retrieve cluster labels/mixing proportions
  if(clust.ind) {
    label.miss   <- missing(labels)
    if(!label.miss)   {
      if(!exists(as.character(substitute(labels)),
          envir=.GlobalEnv))      stop(paste0("Object ", match.call()$labels, " not found"))
      if(length(labels) != n.obs) stop(paste0("'labels' must be a factor of length N=",  n.obs))  
    }
    if(sw["mu.sw"])   {
      mus        <- sims[[G.ind]][[Q.ind]]$mu[,,tmp.store, drop=FALSE]                            
    }
    if(sw["l.sw"])    {
      lmats      <- sims[[G.ind]][[Q.ind]]$load
      if(!is.element(method, c("MFA", "OMFA", "IMFA"))) {
        lmats    <- as.array(lmats)
      }
      lmats      <- lmats[,,,tmp.store, drop=FALSE]
    }
    if(sw["psi.sw"])  {
      psis       <- sims[[G.ind]][[Q.ind]]$psi[,,tmp.store, drop=FALSE]
    }
    if(sw["pi.sw"])   {
      pies       <- sims[[G.ind]][[Q.ind]]$pi.prop[,tmp.store, drop=FALSE]
    }
    z            <- as.matrix(sims[[G.ind]][[Q.ind]]$z.store[,tmp.store])
    if(!label.switch) {
      z.temp     <- factor(z[,1], labels=Gseq)
      for(sl in seq_along(tmp.store)[-1]) {
        sw.lab   <- lab.switch(z.new=z[,sl], z.old=z.temp, Gs=Gseq)
        z[,sl]   <- sw.lab$z
        z.perm   <- sw.lab$z.perm
        if(!identical(as.integer(z.perm), Gseq)) {
          if(sw["mu.sw"])  {
            mus[,Gseq,sl]      <- mus[,z.perm,sl]
          }
          if(sw["l.sw"])   {
            lmats[,,Gseq,sl]   <- lmats[,,z.perm,sl]
          }
          if(sw["psi.sw"]) {
            psis[,Gseq,sl]     <- psis[,z.perm,sl]
          }
          if(sw["pi.sw"])  {
            pies[Gseq,sl]      <- pies[z.perm,sl]
          }
          if(inf.Q)        {
            Q.store[Gseq,sl]   <- Q.store[z.perm,sl]
          }  
        }
      }
    }
    post.z       <- droplevels(setNames(apply(z, 1, function(x) factor(which.max(tabulate(x)), levels=Gseq)), seq_len(n.obs)))
    uncertain    <- 1 - apply(matrix(apply(z, 1, tabulate, nbins=G)/length(tmp.store), nr=nlevels(post.z), nc=n.obs), 2, max)
    if(sw["pi.sw"])    {
      pi.prop    <- pies[Gseq,seq_along(tmp.store), drop=FALSE]
      var.pi     <- apply(pi.prop, 1, var)
      ci.pi      <- apply(pi.prop, 1, function(x) quantile(x, conf.levels))
      post.pi    <- rowMeans(pi.prop, dims=1)
    } else {
      post.pi    <- setNames(prop.table(tabulate(post.z, nbins=G)), paste0("Group ", Gseq))
    }
    if(!label.miss) {
      zlabels    <- factor(labels, labels=seq_along(unique(labels)))
      levs       <- levels(zlabels)
      if(nlevels(post.z) == length(levs)) {
        sw.lab   <- lab.switch(z.new=post.z, z.old=zlabels, Gs=Gseq)
        post.z   <- factor(sw.lab$z)
        l.perm   <- sw.lab$z.perm
        if(sw["mu.sw"])    mus <- mus[,l.perm,, drop=FALSE]
        if(sw["l.sw"])   lmats <- lmats[,,l.perm,, drop=FALSE]
        if(sw["psi.sw"])  psis <- psis[,l.perm,, drop=FALSE]
        gnames   <- paste0("Group ",  l.perm)
        index    <- order(gnames)
        post.pi  <- setNames(post.pi[index], gnames[index])
        if(sw["pi.sw"]) {  
         rownames(pi.prop)     <- gnames
         pi.prop <- pi.prop[index,, drop=FALSE]
         var.pi  <- setNames(var.pi[index],  gnames[index])
         colnames(ci.pi)       <- gnames
         ci.pi   <- ci.pi[,index,   drop=FALSE]
        }
        if(inf.Q)   {
         rownames(Q.store)     <- gnames
         Q.store <- Q.store[index,, drop=FALSE]
        }
      }
      tab        <- table(post.z, zlabels, dnn=list("Predicted", "Observed"))
      tab.stat   <- c(classAgreement(tab), classError(post.z, zlabels))
      if(nrow(tab) != ncol(tab))     {
        tab.stat <- tab.stat[-seq_len(2)]
        names(tab.stat)[4]     <- "error.rate"
      } else {
        names(tab.stat)[6]     <- "error.rate"
      }
      if(tab.stat$error.rate   == 0) {
        tab.stat$misclassified <- NULL
      }
      tab.stat   <- c(list(confusion.matrix = tab), tab.stat)
      uncert.obs <- which(uncertain >= 1/G)
      attr(uncertain, "Obs")   <- if(sum(uncert.obs) != 0) uncert.obs
      if(!label.miss && (nlevels(post.z) == length(levs))) {
        names(tab.stat)[1]     <- "matched.confusion.matrix"
      }
      class(tab.stat)          <- "listof"
    }
    if(isTRUE(MH.step)) {
      alpha.pi   <- sims[[G.ind]][[Q.ind]]$alpha
      post.alpha <- mean(alpha.pi)
      var.alpha  <- var(alpha.pi)
      ci.alpha   <- quantile(alpha.pi, conf.levels)
      rate       <- sims[[G.ind]][[Q.ind]]$rate
      MH.alpha   <- list(alpha.pi = alpha.pi, post.alpha = post.alpha, var.alpha = var.alpha, ci.alpha = ci.alpha, acceptance.rate = rate)
      class(MH.alpha)          <- "listof"
    }
    cluster      <- list(clustering = post.z, z = z, uncertainty = uncertain)
    cluster      <- c(cluster, list(post.pi = post.pi/sum(post.pi)), if(sw["pi.sw"]) list(pi.prop = pi.prop, var.pi = var.pi, 
                      ci.pi = ci.pi), if(!label.miss) list(perf = tab.stat), if(isTRUE(MH.step)) list(MH.alpha = MH.alpha))
    attr(cluster, "Z.init")    <- attr(sims[[G.ind]], "Z.init")
    attr(cluster, "Init.Meth") <- attr(sims, "Init.Z")
    attr(cluster, "Label.Sup") <- !label.miss
    post.z       <- as.numeric(levels(post.z))[post.z]
    z.ind        <- lapply(Gseq, function(g) post.z == g)
  }
  if(inf.Q)   {
    G1           <- G > 1
    Q.tab        <- if(G1) lapply(apply(Q.store, 1, function(x) list(table(x, dnn=NULL))), "[[", 1) else table(Q.store, dnn=NULL)
    Q.prob       <- if(G1) lapply(Q.tab, prop.table) else prop.table(Q.tab)
    Q.mode       <- if(G1) unlist(lapply(Q.tab, function(qt) as.numeric(names(qt[qt == max(qt)])[1]))) else as.numeric(names(Q.tab[Q.tab == max(Q.tab)])[1])
    Q.med        <- if(G1) ceiling(apply(Q.store, 1, median) * 2)/2 else ceiling(median(Q.store) * 2)/2
    if(!Q.T)  {
      Q          <- if(Q.meth == "Mode") Q.mode else floor(Q.med)
    } else    {
      Q          <- if(G.T) Q else setNames(rep(Q, G), paste0("Group ", Gseq))
    }
    if(any(unlist(Q) >= n.var))   warning(paste0("Estimate of Q is not less than the number of variables ", n.var, ": solution may be invalid"), call.=FALSE)
    Q.CI         <- if(G1) apply(Q.store, 1, function(qs) round(quantile(qs, conf.levels))) else round(quantile(Q.store, conf.levels))
    GQ.temp4     <- list(Q = Q, Q.Mode = Q.mode, Q.Median = Q.med, 
                         Q.CI = Q.CI, Q.Probs = Q.prob, Q.Counts = Q.tab)
    GQ.res       <- if(inf.G) c(GQ.temp1, GQ.temp4) else c(list(G = G), GQ.temp4)
    GQ.res       <- c(GQ.res, GQ.temp2)
  }

# Retrieve (unrotated) scores
  no.score       <- all(Q == 0)
  if(no.score)   { 
    if(sw["f.sw"])                warning("Scores & loadings not stored as model has zero factors", call.=FALSE)
    sw["f.sw"]   <- FALSE
  }
  if(sw["f.sw"]) {
    Q.max        <- max(Q) 
    Q.maxs       <- seq_len(Q.max)
    f            <- sims[[G.ind]][[Q.ind]]$f
    if(inf.Q) {
      f          <- as.array(f)
    }
    f            <- f[,Q.maxs,tmp.store, drop=FALSE]
  }

# Loop over g in G to extract other results
  result         <- list(list())
  f.store        <- list()
  mse   <- rmse  <- nrmse  <- cvrmse  <- 
  mad   <- emp.T <- est.T  <- rep(NA, G)
  for(g in Gseq) {
    Qg           <- Q[g]
    Qgs          <- seq_len(Qg)
    sw["l.sw"]   <- attr(sims, "Switch")["l.sw"]
    if(Qg == 0)  {
      if(all(sw["l.sw"],
             !no.score))          warning(paste0("Loadings ", ifelse(G > 1, paste0("for group ", g, " not stored as it"), " not stored as model"), " has zero factors"), call.=FALSE)
      sw["l.sw"] <- FALSE
    }
    store        <- if(inf.Q) seq_along(tmp.store)[which(Q.store[g,] >= Qg)] else seq_along(tmp.store)
    n.store      <- length(store)
  
  # Retrieve (unrotated) loadings  
    if(sw["l.sw"]) {
      if(clust.ind)  {
        lmat     <- adrop(lmats[,Qgs,g,store, drop=FALSE], drop=3)
        l.temp   <- adrop(lmat[,,1, drop=FALSE], drop=3)
      } else {
        lmat     <- sims[[G.ind]][[Q.ind]]$load
        if(inf.Q) {
          lmat   <- as.array(lmat)
        }
        lmat     <- lmat[,Qgs,store, drop=FALSE]
        l.temp   <- adrop(lmat[,,1, drop=FALSE], drop=3)
      }
    }
  
  # Loadings matrix / identifiability / error metrics / etc.  
    if(all(sw["f.sw"], clust.ind)) {
      fg         <- f[z.ind[[g]],Qgs,, drop=FALSE]
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
        if(!data.x)               warning(paste0("Object ", dat, " not found in .GlobalEnv: can't compute empirical covariance and error metrics"), call.=FALSE) 
        dat      <- as.data.frame(get(dat))
        dat      <- dat[sapply(dat, is.numeric)]
        dat      <- scale(dat, center=cent, scale=scaling)
        varnames <- colnames(dat)
      }
      cov.emp    <- cov(dat[z.ind[[g]],, drop=FALSE])
      dimnames(cov.emp)  <- list(varnames, varnames)
      if(sum(z.ind[[g]]) <= 1)    rm(cov.emp)
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
    emp.T[g]     <- exists("cov.emp", envir=environment())
  
  # Compute posterior means and % variation explained
    if(sw["mu.sw"])  {
      post.mu    <- rowMeans(mu, dims=1)
      var.mu     <- apply(mu, 1, var)
      ci.mu      <- apply(mu, 1, function(x) quantile(x, conf.levels))
    }
    if(sw["psi.sw"]) {
      post.psi   <- rowMeans(psi, dims=1)
      var.psi    <- apply(psi, 1, var)
      ci.psi     <- apply(psi, 1, function(x) quantile(x, conf.levels))
    }
    if(sw["l.sw"])   { 
      post.load  <- rowMeans(lmat, dims=2)
      var.load   <- apply(lmat, c(1, 2), var)
      ci.load    <- apply(lmat, c(1, 2), function(x) quantile(x, conf.levels))
      var.exp    <- sum(colSums(post.load * post.load))/n.var
      class(post.load)   <- "loadings"
    } else if(emp.T[g]) {
      var.exp    <- ifelse(sum(z.ind[[g]]) == 0, 0, max(0, (sum(diag(cov.emp)) - sum(post.psi))/n.var))
    }
  
  # Calculate estimated covariance matrices & compute error metrics
    if(clust.ind) {
      if(all(sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        cov.est  <- if(Qg > 0) tcrossprod(post.load) + diag(post.psi) else diag(post.psi)
        if(data.x)      {
          dimnames(cov.est)    <- list(varnames, varnames)
        }
      } else if(g == 1) {
        if(all(!sw["l.sw"], Qg  > 0, !sw["psi.sw"]))  {
                                  warning("Loadings & Uniquenesses not stored: can't estimate covariance matrix and compute error metrics", call.=FALSE)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't estimate covariance matrix and compute error metrics", call.=FALSE)
        } else if(!sw["psi.sw"])  warning("Uniquenesses not stored: can't estimate covariance matrix and compute error metrics", call.=FALSE)
      }  
    } else     {
      cov.est    <- sims[[G.ind]][[Q.ind]]$cov.est
      if(all(recomp, sw["psi.sw"], any(sw["l.sw"], Qg == 0))) {
        cov.est  <- replace(cov.est, is.numeric(cov.est), 0)
        for(r in seq_len(n.store))    {
         sigma   <- if(Qg > 0) tcrossprod(lmat[,,r]) + diag(psi[,r]) else diag(psi[,r])
         cov.est <- cov.est + sigma/n.store
        }
      } else if(all(recomp,  g == 1)) {
        if(all(!sw["l.sw"], Qg  > 0, !sw["psi.sw"]))  {
                                  warning("Loadings & Uniquenesses not stored: can't re-estimate covariance matrix", call.=FALSE)
        } else if(all(Qg > 0,
                  !sw["l.sw"])) { warning("Loadings not stored: can't re-estimate covariance matrix", call.=FALSE)
        } else if(!sw["psi.sw"])  warning("Uniquenesses not stored: can't re-estimate covariance matrix", call.=FALSE) 
      }
    }
    est.T[g]     <- exists("cov.est", envir=environment())
    
    if(all(emp.T[g], est.T[g])) {
      error      <- cov.emp - cov.est
      mse[g]     <- mean(error * error)
      rmse[g]    <- sqrt(mse[g])
      nrmse[g]   <- rmse[g]/(max(cov.emp) - min(cov.emp))
      cvrmse[g]  <- rmse[g]/mean(cov.emp)
      mad[g]     <- mean(abs(error))
      if(any(all(scal.meth != "none", cent) && 
                 sum(round(diag(cov.est))   != 
                 round(diag(cov.emp)))      != 0,
         sum(abs(post.psi  - (1 - post.psi)) < 0) != 0,
         var.exp  > 1))           warning(paste0(ifelse(G == 1, "C", paste0("Group ", g, "'s c")), "hain may not have fully converged"), call.=FALSE)
    }
    
    results      <- list(if(sw["mu.sw"])  list(means     = mu,
                                               var.mu    = var.mu,
                                               ci.mu     = ci.mu), 
                         if(sw["l.sw"])   list(loadings  = lmat, 
                                               post.load = post.load,
                                               var.load  = var.load,
                                               ci.load   = ci.load),
                         if(sw["psi.sw"]) list(psi       = psi,
                                               var.psi   = var.psi,
                                               ci.psi    = ci.psi),
                         if(sw.mx)        list(post.mu   = post.mu), 
                         if(sw.px)        list(post.psi  = post.psi),
                         if(any(sw["l.sw"], 
                                sw.px))   list(var.exp   = var.exp),
                         if(emp.T[g])     list(cov.emp   = cov.emp), 
                         if(est.T[g])     list(cov.est   = cov.est))
    result[[g]]  <- unlist(results, recursive=FALSE)
    attr(result[[g]], "Store") <- n.store
    f.store[[g]] <- store
  }
  if(sw["f.sw"])   {
    f            <- f[,,unique(unlist(f.store)), drop=FALSE]
    scores       <- list(f = f, post.f = rowMeans(f, dims=2), var.f = apply(f, c(1, 2), var),
                         ci.f  = apply(f, c(1, 2), function(x) quantile(x, conf.levels)))
  }
  names(result)  <- paste0("Group", Gseq)
  class(GQ.res)                <- "listof"
  attr(GQ.res, "Criterion")    <- criterion
  attr(GQ.res, "Factors")      <- n.fac
  attr(GQ.res, "Groups")       <- n.grp
  attr(GQ.res, "Supplied")     <- c(Q=Q.T, G=G.T)
  err.T                        <- unlist(lapply(Gseq, function(g) all(emp.T[g], est.T[g])))
  if(any(err.T))   {
    errors       <- lapply(list(MSE = mse, RMSE = rmse, NRMSE = nrmse, CVRMSE = cvrmse, MAD = mad), setNames, paste0("Group ", Gseq))
    if(G > 1)      {
      errors     <- c(errors, list(Averages = unlist(lapply(errors, mean, na.rm=TRUE))))
      class(errors)            <- "listof"
    } else {
      errors     <- setNames(unlist(errors), names(errors))
    }
  }
  
  if(sw["mu.sw"])  {
    post.mu      <- do.call(rbind, lapply(result, "[[", "post.mu"))
    var.mu       <- do.call(rbind, lapply(result, "[[", "var.mu"))
    ci.mu        <- Filter(Negate(is.null), lapply(result, "[[", "ci.mu"))  
    means        <- list(post.mu = post.mu, var.mu = var.mu, ci.mu = ci.mu)
  }
  if(sw["l.sw"])   {
    post.load    <- Filter(Negate(is.null), lapply(result, "[[", "post.load"))
    var.load     <- Filter(Negate(is.null), lapply(result, "[[", "var.load"))
    ci.load      <- Filter(Negate(is.null), lapply(result, "[[", "ci.load"))  
    loads        <- list(post.load = post.load, var.load = var.load, ci.load = ci.load)
  }
  if(sw["psi.sw"]) {
    post.psi     <- do.call(rbind, lapply(result, "[[", "post.psi"))
    var.psi      <- do.call(rbind, lapply(result, "[[", "var.psi"))
    ci.psi       <- Filter(Negate(is.null), lapply(result, "[[", "ci.psi"))  
    psis         <- list(post.psi = post.psi, var.psi = var.psi, ci.psi = ci.psi)
  }
  result         <- c(result, if(exists("cluster", envir=environment())) list(Clust = cluster), 
                      if(any(err.T))   list(Error        = errors),  list(GQ.results = GQ.res), 
                      if(sw["mu.sw"])  list(Means        =  means),
                      if(sw["l.sw"])   list(Loadings     =  loads),
                      if(sw["f.sw"])   list(Scores       = scores),
                      if(sw["psi.sw"]) list(Uniquenesses =   psis))
                      
  class(result)                <- "IMIFA"
  attr(result, "Conf.Level")   <- conf.level
  attr(result, "Errors")       <- any(err.T)
  attr(result, "Method")       <- method
  if(is.element(method, c("IMFA", "IMIFA"))) {
    attr(result, "MH.step")    <- MH.step
    attr(result, "Gen.Slice")  <- attr(sims, "Gen.Slice")
  }
  attr(result, "Name")         <- attr(sims, "Name")
  attr(result, "Obs")          <- n.obs
  attr(result, "Store")        <- tmp.store
  attr(result, "Switch")       <- sw
  attr(result, "Vars")         <- n.var
  return(result)
}