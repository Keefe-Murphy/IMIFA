#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages    <- c("abind", "e1071", "matrixStats", "mclust", "MCMCpack", "plotrix", "slam")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  suppressMessages(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  suppressMessages(lapply(setdiff(packages, (.packages())), library, ch=T))
}
rm(packages)
message("   ________  __________________\n  /_  __/  |/   /_  __/ ___/ _ \\  \n   / / / /|_// / / / / /__/ /_\\ \\ \n _/ /_/ /   / /_/ /_/ ___/ /___\\ \\ \n/____/_/   /_/_____/_/  /_/     \\_\\    version 1.0")

imifa.mcmc  <- function(dat = NULL, method = c("IMIFA", "MIFA", "MFA", "IFA", "FA", "classify"), 
                        n.iters = 50000, Labels = NULL, factanal = F, Q.star = NULL, range.G = NULL, 
                        range.Q = NULL, Q.fac = NULL,  burnin = n.iters/5, thinning = 2, centering = T, 
                        scaling = c("unit", "pareto", "none"), verbose = F, adapt = T, b0 = NULL, b1 = NULL, 
                        prop = NULL, epsilon = NULL, sigma.mu = NULL, sigma.l = NULL, mu0g = F, psi0g = F, mu.zero = NULL,
                        phi.nu = NULL, psi.alpha = NULL, psi.beta = NULL, alpha.d1 = NULL, alpha.dk = NULL, beta.d1 = NULL,
                        beta.dk = NULL, alpha.pi = NULL, z.init = c("kmeans", "list", "mclust", "priors"), z.list = NULL, 
                        profile = F, mu.switch = T, f.switch = T, load.switch = T, psi.switch = T, pi.switch = T) {
  
  defpar    <- suppressWarnings(par(no.readonly = T))
  defop     <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(method == "classification") {
    method  <- "classify"
  }
  method    <- match.arg(method)
  if(is.logical(scaling) && !scaling) {
    scaling <- "none"
  }
  scaling   <- match.arg(scaling)
  if(missing(dat))                  stop("Dataset must be supplied")
  if(!exists(deparse(substitute(dat)),
             envir=.GlobalEnv))     stop(paste0("Object ", match.call()$dat, " not found"))
  if(!is.logical(factanal))         stop("'factanal' must be TRUE or FALSE")
  if(!is.logical(centering))        stop("'centering' must be TRUE or FALSE")
  if(!is.logical(verbose))          stop("'verbose' must be TRUE or FALSE")
  if(!is.logical(profile))          stop("'profile' must be TRUE or FALSE")
  
# Remove non-numeric columns & apply centering & scaling if necessary 
  burnin    <- as.integer(burnin)
  thinning  <- as.integer(thinning)
  n.iters   <- max(burnin + 1, as.integer(n.iters))
  iters     <- seq(from=burnin + 1, to=n.iters, by=thinning)
  iters     <- iters[iters > 0]
  dat       <- as.data.frame(dat)
  raw.dat   <- dat[sapply(dat, is.numeric)]
  if(scaling != "none") {
    scal    <- apply(raw.dat, 2, sd)
    if(scaling == "pareto") {
      scal  <- sqrt(scal)
    }
  } else {
    scal    <- F
  }
  dat       <- scale(raw.dat, center=centering, scale=scal)
  N         <- nrow(dat)
  P         <- ncol(dat)
  
# Manage storage switches & warnings for other function inputs
  if(!missing(mu.switch) && all(!mu.switch, !centering)) {
                                    warning("Centering hasn't been applied - are you sure you want mu.switch=F?", call.=F)
  }
  switches  <- c(mu.sw=mu.switch, f.sw=f.switch, l.sw=load.switch, psi.sw=psi.switch, pi.sw=pi.switch)
  if(!is.logical(switches))         stop("All logical switches must be TRUE or FALSE")
  if(!is.element(method, c("MFA", "MIFA"))) {
    if(!missing(range.G) &&  
       any(range.G  > 1))           warning(paste0("'range.G' must be 1 for the ", method, " method"), call.=F)
    range.G <- 1
    meth    <- method
  } else {
    if(missing(range.G))            stop("'range.G' must be specified")
    if(any(range.G  < 1))           stop("'range.G' must be strictly positive")
    range.G <- sort(unique(range.G))
    meth    <- rep(method, length(range.G))
    if(range.G[1]  == 1)  {
      if(meth[1] != "IMIFA")  {                                    
        if(meth[1] == "MFA")  {
          meth[1]  <- "FA"                              
        } 
        if(meth[1] == "MIFA") {
          meth[1]  <- "IFA"
        } 
      }
                                    warning(paste0("Forced use of ", meth[1], " method where range.G is equal to 1"), call.=F)
    }                               
  }
  no.fac    <- is.element(method, c("FA", "MFA")) && all(range.Q == 0)
  if(is.element(method, c("FA", "MFA"))) {
    if(!missing(Q.star))  {
      rm(Q.star)
                                    warning(paste0("'Q.star' is not used for the ", method, " method"), call.=F)
    }            
    if(missing(range.Q))            stop("'range.Q' must be specified")
    if(any(range.Q < 0))            stop("'range.Q' must be strictly non-negative")
    range.Q <- sort(unique(range.Q))   
  } else {
    if(!missing(range.Q)) {
      rm(range.Q)
                                    warning(paste0("'range.Q' is not used for the ", method, " method"), call.=F)
    }           
    if(missing(Q.star)) {
      Q.star       <- min(floor(3 * log(P)), P, N - 1)
    } else {
      if(Q.star    <= 0)            stop("Q.star must be strictly positive")
      if(Q.star     > P)            stop(paste0("Number of factors must be less than the number of variables, ", P))
      if(Q.star    >= N)            stop(paste0("Number of factors must be less than the number of observations, ", N))
    } 
    if(!is.logical(adapt))          stop("'adapt' must be TRUE or FALSE") 
  }
  if(no.fac) {   
    if(all(switches[c("f.sw", "l.sw")])) {
                                    warning("Scores & Loadings not stored as model has zero factors", call.=F)
    } else if(switches["f.sw"])   { warning("Scores not stored as model has zero factors", call.=F)
    } else if(switches["l.sw"])   { warning("Loadings not stored as model has zero factors", call.=F)
    }                               
    switches[c("f.sw", "l.sw")]  <- F                              
  } else {
    if(all(!switches[c("f.sw", "l.sw")])) { 
                                    warning("Posterior Scores & Loadings won't be available as they're not being stored", call.=F)
    } else if(!switches["f.sw"])  { warning("Posterior Scores won't be available as they're not being stored", call.=F)
    } else if(!switches["l.sw"])  { warning("Posterior Loadings won't be available as they're not being stored", call.=F)
    }
  }

  if(any(all(method == "MFA",  any(range.G > 1)) && any(range.Q > 0),
         all(method == "MIFA", any(range.G > 1)),
             method == "IMIFA"))  {
    if(all(!switches["l.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Loadings & Psi not stored: unable to estimate covariance matrix and compute error metrics", call.=F)
    } else if(!switches["l.sw"])  { warning("Loadings not stored: unable to estimate covariance matrix and compute error metrics", call.=F)
    } else if(!switches["psi.sw"])  warning("Psi not stored: unable to estimate covariance matrix and compute error metrics", call.=F)
  }
  if(any(all(method == "MFA",  any(range.G > 1)),
         all(method == "MIFA", any(range.G > 1)),
         method == "IMIFA"))      {
    if(all(!switches["mu.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Means & Psi not stored: posterior mean estimates won't be available", call.=F)
    } else if(!switches["mu.sw"]) { warning("Means not stored: posterior mean estimates won't be available", call.=F)
    } else if(!switches["psi.sw"])  warning("Psi not stored: posterior mean estimates won't be available", call.=F)
  }
  
# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  cov.mat          <- var(dat)
  if(is.null(rownames(dat))) rownames(dat) <- seq_len(N)
  if(missing("sigma.mu"))    sigma.mu      <- diag(cov.mat)
  if(scaling == "unit")      sigma.mu      <- sigma.mu[1]
  if(any(sigma.mu <= 0))            stop("'sigma.mu' must be strictly positive")
  if(missing("psi.alpha"))   psi.alpha     <- 2.5
  if(psi.alpha <= 0)                stop("'psi.alpha' must be strictly positive")
  if(is.element(method, c("FA", "MFA"))) {
    if(missing("sigma.l"))   sigma.l       <- 0.5
    if(sigma.l <= 0)                stop("'sigma.l' must be strictly positive")
  } else {
    if(missing("phi.nu"))    phi.nu        <- 1.5
    if(phi.nu <= 0)                 stop("'phi.nu' must be strictly positive")
    if(missing("alpha.d1"))  alpha.d1      <- 2
    if(alpha.d1 <= 0)               stop("'alpha.d1' must be strictly positive")
    if(missing("alpha.dk"))  alpha.dk      <- 10
    if(alpha.dk <= 1)               stop("'alpha.dk' must be greater than 1")
    if(missing("beta.d1"))   beta.d1       <- 1
    if(beta.d1 <= 0)                stop("'beta.d1' must be strictly positive")
    if(missing("beta.dk"))   beta.dk       <- 1
    if(beta.dk <= 0)                stop("'beta.dk' must be strictly positive")
    if(missing("b0"))        b0            <- 0.1
    if(b0  < 0)                     stop("'b0' must be positive to ensure valid adaptation probability")
    if(missing("b1"))        b1            <- 0.00005
    if(b1 <= 0)                     stop("'b1' must be strictly positive to ensure adaptation probability decreases")
    if(missing("prop"))      prop          <- 3/4
    if(abs(prop - (1 - prop)) < 0)  stop("'prop' must be a single number between 0 and 1")
    if(missing("epsilon"))   epsilon       <- ifelse(centering, 0.1, 0.005)
    if(abs(epsilon - 
          (1 - epsilon)) < 0)       stop("'epsilon' must be a single number between 0 and 1")
  } 
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(!is.logical(mu0g))           stop("'mu0g' must be TRUE or FALSE")
    if(!is.logical(psi0g))          stop("'psi0g' must be TRUE or FALSE")
    if(missing("alpha.pi"))  alpha.pi      <- ifelse(method == "IMIFA", 0.1, 0.5)
    if(abs(alpha.pi -
          (1 - alpha.pi)) < 0)      stop("'alpha.pi' must be a single number between 0 and 1")
                             z.init        <- match.arg(z.init)
    if(method != "IMIFA") {
      if(!missing(z.list))   {
        if(!is.list(z.list))   z.list      <- list(z.list)
                               z.list      <- lapply(z.list, as.factor)
        if(z.init != "list") { z.init      <- "list"
                                    warning("'z.init' set to 'list' as 'z.list' was supplied", call.=F) }
        if(length(z.list)   != length(range.G))      {
                                    stop(paste0("'z.list' must be a list of length ", length(range.G))) }
        if(!all(lapply(z.list, nlevels) == range.G)) {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'")) }
        if(!all(lapply(z.list, length)  == N))       {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N)) }
      }
      if(all(missing(z.list),  z.init   == "list"))  {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'")) }
    }
  }
  if(all(!is.element(method, c("MFA", "MIFA")), 
         !missing(z.init) || 
         !missing(z.list)))         warning(paste0("z does not need to be initialised for the ", method, " method"), call.=F)
  
  if(method == "classify") {
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", "IFA", ".R", sep=""), local=T)
  } else {
    for(g.meth in unique(meth)) {
      source(paste(getwd(), "/IMIFA-GIT/Gibbs_", g.meth, ".R", sep=""), local=T)
    }
  }
  source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=T)

  imifa     <- list(list())
  Gi        <- 1
  Qi        <- 1
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, 
                    thinning = thinning, iters = iters, verbose = verbose, sw = switches)
  if(!is.element(method, c("FA", "MFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(phi.nu = phi.nu, alpha.d1 = alpha.d1, alpha.dk = alpha.dk, beta.d1 = beta.d1,
                                             beta.dk = beta.dk, adapt = adapt, b0 = b0, b1 = b1, prop = prop, epsilon = epsilon))
  } else {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l = sigma.l))
  }

  init.start       <- proc.time()  
  beta.x           <- missing("psi.beta")
  mu0.x            <- missing("mu.zero")
  mu               <- list(colMeans(dat))
  if(beta.x) {
    psi.beta       <- temp.psi <- list(psi.hyper(psi.alpha, cov.mat))
  } else {
    psi.beta       <- len.check(psi.beta)
  }
  if(mu0.x)  {
    mu.zero        <- mu
  } else {
    mu.zero        <- len.check(mu.zero)
  }
  if(is.element(method, c("MFA", "MIFA"))) {
    clust          <- list()
    pi.alpha       <- list()
    pi.prop        <- list()
    zi             <- list()
    for(g in seq_along(range.G)) {
      G            <- range.G[g]
      pi.alpha[[g]]     <- rep(alpha.pi, G)
      if(z.init    == "kmeans")     {
        k.res      <- try(kmeans(dat, G, nstart=100), silent=T)
        if(!inherits(k.res, "try-error"))  {
          zi[[g]]  <- as.numeric(factor(k.res$cluster, levels=seq_len(G)))
        } else                      warning("Cannot initialise cluster labels using kmeans. Try another z.init method")
      } else if(z.init  == "list")   {
        zi[[g]]    <- as.numeric(z.list[[g]])
      } else if(z.init  == "mclust") {
        m.res      <- try(Mclust(dat, G), silent=T)
        if(!inherits(m.res, "try_error"))  {
          zi[[g]]  <- as.numeric(m.res$classification)
        } else                      warning("Cannot initialise cluster labels using mclust. Try another z.init method")
      } else {
        zips       <- rep(1, N)
        while(all(length(unique(zips)) != G,
              any(prop.table(tabulate(zips, nbins=G)) < 1/G^2))) {
          pi.prop  <- sim.pi(pi.alpha=pi.alpha[[g]], nn=0)
          zips     <- sim.z.p(N=N, prob.z=pi.prop)
        }
        zi[[g]]    <- as.numeric(zips)
        rm(zips)
      }
      pi.prop[[g]] <- t(prop.table(tabulate(zi[[g]], nbins=G)))
      mu[[g]]      <- do.call(cbind, lapply(seq_len(G), function(gg) if(pi.prop[[g]][,gg] > 0) colMeans(dat[zi[[g]] == gg,, drop=F]) else rep(0, P)))
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]] else do.call(cbind, lapply(seq_len(G), function(gg) colMeans(dat)))  
      }
      if(beta.x)  {
        if(psi0g) {
          cov.gg   <- lapply(seq_len(G), function(gg) if(pi.prop[[g]][,gg] > 0) cov(dat[zi[[g]] == gg,, drop=F]) else cov.mat)
          psi.beta[[g]] <- do.call(cbind, lapply(seq_len(G), function(gg) psi.hyper(psi.alpha, cov.gg[[gg]])))
        } else {
          psi.beta[[g]] <- replicate(G, temp.psi[[1]])
        }
      }
      clust[[g]]   <- list(z = zi[[g]], pi.alpha = pi.alpha[[g]], 
                           pi.prop = pi.prop[[g]], label.switch = c(mu0g, psi0g))
    }
  }
  if(all(round(unlist(sapply(mu.zero, sum))) == 0)) {
    mu.zero        <- lapply(mu.zero, function(x) 0)
  }
  if(any(unlist(sapply(psi.beta, 
         function(x) sapply(x, 
         function(x) x <= 0)))))    stop("'psi.beta' must be strictly positive")
  init.time        <- proc.time() - init.start
  
  if(profile)  Rprof()
  if(is.element(method, c("IFA", "MIFA"))) {
    if(length(range.G) == 1) {
      start.time   <- proc.time()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),                          
                                     args=append(list(data = dat, N = N, G = range.G, Q = Q.star, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MIFA") clust[[Gi]]), gibbs.arg))
    } else {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = Q.star, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MIFA") clust[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Gi/length(range.G) * 100, 2), "% Complete\n"))
      }
    }
  } else if(is.element(method, c("FA", "MFA")))   {
    if(any(range.Q >= P))           stop(paste0("Number of factors must be less than the number of variables ", P))
    if(any(range.Q >= N - 1))       stop(paste0("Number of factors must be less than the number of observations minus 1 ", N - 1))
    if(all(length(range.G) == 1, length(range.Q) == 1)) {
      start.time   <- proc.time()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]), 
                                     args=append(list(data = dat, N = N, G = range.G, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MFA") clust[[Gi]]), gibbs.arg))
    } else if(length(range.G) == 1) {
      start.time   <- proc.time()
      for(q in range.Q) { 
        Qi         <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = range.G, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MFA") clust[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Qi/length(range.Q) * 100, 2), "% Complete\n"))
      }
    } else if(length(range.Q) == 1) {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MFA") clust[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Gi/length(range.G) * 100, 2), "% Complete\n"))
      }
    } else {
      counter      <- 0
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        for(q in range.Q) {
          Qi       <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MFA") clust[[Gi]]), gibbs.arg))
        counter    <- counter + 1
        if(verbose)                 cat(paste0(round(counter/(length(range.G) * length(range.Q)) * 100, 2), "% Complete\n"))
        }
      }
    }
  } else if(method == "classify") {
    if(missing(Labels))             stop("Data must be labelled for classification")
    if(!exists(deparse(substitute(Labels)),
               envir=.GlobalEnv))   stop(paste0("Object ", match.call()$Labels, " not found"))
    Labels  <- as.factor(Labels)
    if(length(Labels) != N)         stop(paste0("Labels must be a factor of length N=",  n.obs))
    range.G        <- nlevels(Labels)
    start.time     <- proc.time()
    for(g in seq_len(range.G)) {
      temp.dat     <- dat[Labels == levels(Labels)[g],]
      imifa[[g]]          <- list()
      imifa[[g]][[Qi]]    <- do.call(paste0("gibbs.", "IFA"),
                                     args=append(list(data = temp.dat, N = nrow(temp.dat), Q = Q.star), gibbs.arg))
      if(verbose)                   cat(paste0(round(g/range.G * 100, 2), "% Complete\n"))
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "MFA", length(range.G) * length(range.Q),
                          ifelse(method == "FA",  length(range.Q), 
                            ifelse(method == "MIFA", length(range.G),
                              ifelse(method == "classify", nlevels(Labels), 
                                     length(Q.star)))))
  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }
  dat.name  <- as.character(match.call()$dat)
  if(is.element(method, c("FA", "MFA")))   {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  } else {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], "IFA"))
  }
  if(is.element(method, c("MFA", "MIFA"))) {
    for(g in seq_along(range.G)) {
      attr(imifa[[g]], 
           "Z.init")      <- factor(zi[[g]], levels=seq_len(range.G[g]))
    }
  }
  gnames    <- paste0(range.G, ifelse(range.G == 1, "Group", "Groups"))
  names(imifa)            <- gnames
  attr(imifa, "Center")   <- centering
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors")  <- if(is.element(method, c("FA", "MFA"))) range.Q else Q.star
  attr(imifa, "Groups")   <- if(method != "IMIFA") range.G else G.star
  if(is.element(method, c("MFA", "MIFA"))) {
    attr(imifa, "Init.Z") <- z.init
    attr(imifa, 
         "Label.Switch")  <- any(mu0g, psi0g)
  }
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.name
  attr(imifa, "Obs")      <- N
  attr(imifa, "Scaling")  <- scal
  attr(attr(imifa,
  "Scaling"), "Method")   <- scaling
  attr(imifa, "Store")    <- length(iters)
  attr(imifa, "Switch")   <- switches
  if(any(is.element(method, c("IFA", "IMIFA")), 
     (is.element(method, c("FA", "MFA")) && length(range.Q) == 1)))       {
    attr(imifa, "Time")   <- round(tot.time, 2)
  } else if(all(is.element(method, c("MFA", "MIFA")), z.init  != "list")) {
    attr(imifa, "Time")   <- lapply(list(Total = tot.time, Average = avg.time, Z.Initialisation = init.time), function(x) round(x, 2)) 
  } else {
    attr(imifa, "Time")   <- lapply(list(Total = tot.time, Average = avg.time), function(x) round(x, 2)) 
  }
  attr(imifa, "Vars")     <- P
  if(verbose)                print(attr(imifa, "Time"))  
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac))    factanal <- T
  if(factanal) {
    if(missing(Q.fac)) {
      if(missing(range.Q)) {
        Q.fac      <- round(sqrt(P))
      } else {
        Q.fac      <- max(1, max(range.Q))
      }
    }
    fac     <- try(factanal(dat, factors=Q.fac, control=list(nstart=50)))
    if(!inherits(fac, "try-error")) {
      imifa <- append(imifa, list(fac = fac))
      names(imifa)[length(imifa)] <- "Factanal"
    }
  } 
  class(imifa)     <- "IMIFA"
  return(imifa)
}

source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/SimulateData.R", sep=""))