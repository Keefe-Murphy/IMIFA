#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages    <- c("abind", "dichromat", "e1071", "gclus", "matrixStats", "mclust", "MCMCpack", "mvnfast", "plotrix", "slam")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  suppressMessages(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  suppressMessages(lapply(setdiff(packages, (.packages())), library, ch=TRUE))
}
rm(packages)
if(!exists("mcmc.IMIFA", 
    envir=.GlobalEnv))              packageStartupMessage("   ________  __________________\n  /_  __/  |/   /_  __/ ___/ _ \\  \n   / / / /|_// / / / / /__/ /_\\ \\ \n _/ /_/ /   / /_/ /_/ ___/ /___\\ \\ \n/____/_/   /_/_____/_/  /_/     \\_\\    version 1.0")
source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/SimulateData.R", sep=""))

mcmc.IMIFA  <- function(dat = NULL, method = c("IMIFA", "IMFA", "OMIFA", "OMFA", "MIFA", "MFA", "IFA", "FA", "classify"), 
                        n.iters = 50000, Labels = NULL, factanal = FALSE, range.G = NULL, range.Q = NULL, verbose = FALSE, Q.fac = NULL,  
                        burnin = n.iters/5, thinning = 2, centering = TRUE, scaling = c("unit", "pareto", "none"), trunc.G = NULL, MH.lower = NULL,
                        adapt = TRUE, b0 = NULL, b1 = NULL, delta0g = FALSE, prop = NULL, epsilon = NULL, sigma.mu = NULL, sigma.l = NULL, MH.step = TRUE,
                        mu0g = FALSE, psi0g = FALSE, mu.zero = NULL, phi.nu = NULL, psi.alpha = NULL, psi.beta = NULL, alpha.d1 = NULL, pp = NULL, MH.upper = NULL,
                        alpha.dk = NULL, beta.d1 = NULL, beta.dk = NULL, alpha.pi = NULL, z.list = NULL, profile = FALSE, mu.switch = TRUE, gen.slice = FALSE,
                        f.switch = TRUE, load.switch = TRUE, psi.switch = TRUE, pi.switch = TRUE, z.init = c("kmeans", "list", "mclust", "priors")) {
  
  defpar    <- suppressWarnings(par(no.readonly=TRUE))
  defopt    <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defopt)), add=TRUE)
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
                                    warning("Centering hasn't been applied - are you sure you want mu.switch=FALSE?", call.=FALSE)
  }
  switches  <- c(mu.sw=mu.switch, f.sw=f.switch, l.sw=load.switch, psi.sw=psi.switch, pi.sw=pi.switch)
  if(!is.logical(switches))         stop("All logical switches must be TRUE or FALSE")
  G.x       <- missing(range.G)
  if(!is.element(method, c("MFA", "MIFA")))    {
    if(all(!G.x, is.element(method, c("FA", "IFA"))) &&  
       any(range.G  > 1))           warning(paste0("'range.G' must be 1 for the ", method, " method"), call.=FALSE)
    if(is.element(method, c("OMIFA", "OMFA", "IMFA", "IMIFA"))) {
      if(G.x) {
        range.G    <- max(2, floor(2 * log(N)))
      }
      if(range.G   == 1)            stop(paste0("'range.G' should be at least greater than 1 for the ", method, " method"))
      if(is.element(method, c("IMFA", "IMIFA"))) {
        if(missing(trunc.G))  {
          trunc.G  <- ifelse(N < 100, N, 100)
        } 
        if(!is.logical(gen.slice))  stop("'gen.slice' must be TRUE or FALSE") 
        if(!is.logical(MH.step))    stop("'MH.step' must be TRUE or FALSE") 
        if(missing(pp)) {
          pp       <- 0.5
        }
        if(missing(MH.lower)) {
          MH.lower <- ifelse(N >= P, 0.5, 0)
        }
        if(missing(MH.upper)) {
          MH.upper <- range.G/2
        }
        if(all(length(pp)  > 1,
           pp < 0  && pp   > 1))    stop("'pp' must be a single number between 0 and 1")
        if(all(length(MH.lower) > 1,
           MH.lower < 0))           stop("'MH.lower' must be single number, strictly positive")
        if(all(length(MH.upper) > 1,
           MH.upper < 1))           stop("'MH.upper' must be single number, at least 1")
        if(MH.upper <= MH.lower)    stop(paste0("'MH.upper'(=", MH.upper, ") must be greater than 'MH.lower'(=", MH.lower, ")"))
        if(length(trunc.G) > 1)     stop("'trunc.G' must be a single number")
        if(all(N    > 100, 
           trunc.G  < 100))         stop("'trunc.G' must be at least 100")
        if(trunc.G  < range.G)      stop(paste0("'trunc.G' must be at least range.G=", range.G))
        if(trunc.G  > N)            stop(paste0("'trunc.G' cannot be greater than N=", N))
      }
    } else {
      range.G <- 1
    }
    if(length(range.G) != 1)        stop(paste0("Only one range.G value can be specified for the ", method, " method"))
    meth    <- method
  } else {
    if(G.x)                         stop("'range.G' must be specified")
    if(any(range.G  < 1))           stop("'range.G' must be strictly positive")
    range.G <- sort(unique(range.G))
    meth    <- rep(method, length(range.G))                               
  }
  if(any(range.G >= N))             stop(paste0("'range.G' must be less than the number of observations N=", N))
  if(range.G[1]  == 1)  {
    if(is.element(meth[1], c("IMIFA", "IMFA",
       "OMIFA", "OMFA")))   {       stop("'method' must be FA or IFA for a one group model")
    } else {
      if(meth[1] == "MFA")  {
        meth[1]  <- "FA"                              
      } 
      if(meth[1] == "MIFA") {
        meth[1]  <- "IFA"
      } 
    }
    if(!is.element(method, 
                   c("FA", "IFA"))) message(paste0("Forced use of ", meth[1], " method where 'range.G' is equal to 1"))
  }
  
# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  cov.mat          <- var(dat)
  if(is.null(rownames(dat))) rownames(dat) <- seq_len(N)
  if(missing("sigma.mu"))    sigma.mu      <- diag(cov.mat)
  if(scaling == "unit")      sigma.mu      <- sigma.mu[1]
  if(any(sigma.mu  <= 0))           stop("'sigma.mu' must be strictly positive")
  if(missing("psi.alpha"))   psi.alpha     <- 2.5
  if(psi.alpha <= 0)                stop("'psi.alpha' must be strictly positive")
  Q.miss    <- missing(range.Q)
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(Q.miss)                      stop("'range.Q' must be specified") 
    if(any(range.Q < 0))            stop(paste0("'range.Q' must be non-negative for the ", method, " method"))
  } else {
    if(Q.miss)        range.Q    <- min(floor(3 * log(P)), P, N - 1)
    if(range.Q    <= 0)             stop(paste0("'range.Q' must be strictly positive for the ", method, " method"))
  }
  range.Q   <- sort(unique(range.Q))  
  if(any(range.Q   > P))            stop(paste0("Number of factors must be less than the number of variables, ", P))
  if(any(range.Q  >= N))            stop(paste0("Number of factors must be less than the number of observations, ", N))
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(missing("sigma.l"))   sigma.l       <- 1
    if(sigma.l <= 0)                stop("'sigma.l' must be strictly positive")            
  } else {            
    if(!is.logical(adapt))          stop("'adapt' must be TRUE or FALSE") 
    if(missing("phi.nu"))    phi.nu        <- 1.5
    if(phi.nu  <= 0)                stop("'phi.nu' must be strictly positive")
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
  if(any(all(method == "MFA",  any(range.G > 1)) && any(range.Q > 0),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA",
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["l.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Loadings & Psi not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Loadings not stored: may be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Psi not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
  }
  if(any(all(method == "MFA",  any(range.G > 1)),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA", 
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["mu.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Means & Psi not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["mu.sw"]) { warning("Means not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Psi not stored: posterior mean estimates won't be available", call.=FALSE)
  }
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")) && all(range.Q == 0)) {   
    if(all(switches[c("f.sw", "l.sw")]))  {
                                    warning("Scores & Loadings not stored as model has zero factors", call.=FALSE)
    } else if(switches["f.sw"])   { warning("Scores not stored as model has zero factors", call.=FALSE)
    } else if(switches["l.sw"])   { warning("Loadings not stored as model has zero factors", call.=FALSE)
    }                               
    switches[c("f.sw", "l.sw")]  <- F                              
  } else {
    if(all(!switches[c("f.sw", "l.sw")])) { 
                                    warning("Posterior Scores & Loadings won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["f.sw"])  { warning("Posterior Scores won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Posterior Loadings won't be available as they're not being stored", call.=FALSE)
    }
  }
  if(all(is.element(method, c("FA", "IFA")), 
         !missing(z.init) || 
         !missing(z.list)))         message(paste0("z does not need to be initialised for the ", method, " method"))
  if(!is.logical(mu0g))             stop("'mu0g' must be TRUE or FALSE")
  if(!is.logical(psi0g))            stop("'psi0g' must be TRUE or FALSE")
  if(!is.logical(delta0g))          stop("'delta0g' must be TRUE or FALSE")
  sw0gs     <- c(mu0g = mu0g, psi0g = psi0g, delta0g = delta0g)
  if(all(!is.element(method, c("MFA", "MIFA")), 
         any(sw0gs)))               stop(paste0(names(which(sw0gs)), " should be FALSE for the ", method, " method\n"))
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(all(method != "MIFA",
       delta0g))                    stop("'delta0g' and can only be TRUE for the 'MIFA' method")
    alpha.miss  <- missing("alpha.pi")
    if(alpha.miss)           alpha.pi      <- ifelse(is.element(method, c("OMIFA", "OMFA")), 0.5/range.G, ifelse(all(MH.step, 
                                                     is.element(method, c("IMIFA", "IMFA"))), runif(1, MH.lower, MH.upper), 1))
    if(all(is.element(method,  c("IMIFA", "IMFA")),
           alpha.miss,  !MH.step))  warning("'alpha.pi' fixed at 1 rather than simulated from prior as it's not being learned via Metropolis-Hastings updates", call.=FALSE)
    if(length(alpha.pi) != 1)       stop("'alpha.pi' must be specified as a scalar to ensure an exchangeable prior")
    if(alpha.pi <= 0)               stop("'alpha.pi' must be strictly positive")
    if(all(!is.element(method, c("IMFA", "IMIFA")),
           alpha.pi  > 1))          warning("Are you sure alpha.pi should be greater than 1?", call.=FALSE)
                             z.init        <- match.arg(z.init)
    if(all(is.element(method,  c("OMIFA", "OMFA")), !is.element(z.init, 
       c("list", "kmeans"))))       stop("'z.init' must be set to 'list' or 'kmeans' for the OMIFA method to ensure all groups are populated at the initialisation stage")
    if(!missing(z.list))   {
      if(!is.list(z.list))   z.list        <- list(z.list)
                             z.list        <- lapply(z.list, as.factor)
      if(z.init != "list") { z.init        <- "list"
                                    message("'z.init' set to 'list' as 'z.list' was supplied") }
      if(length(z.list)   != length(range.G))        {
                                    stop(paste0("'z.list' must be a list of length ", length(range.G))) }
                             list.levels   <- lapply(z.list, nlevels)
      if(!all(list.levels == range.G))               {
        if(!is.element(method, c("IMIFA", 
                       "IMFA", "OMIFA", "OMFA")))    {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'")) 
        } else                      stop(paste0("Only ", list.levels, " groups are populated according to z.list, but 'range.G' has been set to ", range.G, ".\n  Reset range.G to this value to avoid redunandtly carrying around empty groups"))
      }
      if(!all(lapply(z.list, length)    == N))       {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N)) }
    }
    if(all(missing(z.list),  z.init     == "list"))  {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'")) }
  }
  if(method == "classify") {
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", "IFA", ".R", sep=""), local=TRUE)
  } else {
    for(g.meth in unique(meth)) {
      source(paste(getwd(), "/IMIFA-GIT/Gibbs_", g.meth, ".R", sep=""), local=TRUE)
    }
  }
  source(paste(getwd(), "/IMIFA-GIT/FullConditionals.R", sep=""), local=TRUE)

  imifa     <- list(list())
  Gi        <- 1
  Qi        <- 1
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, 
                    thinning = thinning, iters = iters, verbose = verbose, sw = switches)
  if(is.element(method, c("IMIFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(trunc.G = trunc.G, pp = pp, gen.slice = gen.slice, MH.step = MH.step, MH.lower = MH.lower, MH.upper = MH.upper))
  }
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(phi.nu = phi.nu, beta.d1 = beta.d1, beta.dk = beta.dk, 
                                             adapt = adapt, b0 = b0, b1 = b1, prop = prop, epsilon = epsilon))
    temp.args      <- gibbs.arg
  } else {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l = sigma.l))
  }

  init.start       <- proc.time()  
  mu               <- list(colMeans(dat))
  beta.x           <- missing("psi.beta")
  mu0.x            <- missing("mu.zero")
  ad1.x            <- all(missing("alpha.d1"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")))
  adk.x            <- all(missing("alpha.dk"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")))
  if(all(z.init != "list", any(sw0gs))) {
    if(delta0g)                     stop("'delta0g' can only be TRUE if z.init=list\n")
    if(all(!mu0.x, mu0g))           stop("'mu.zero' can only be supplied for each group if z.init=list")
    if(all(!beta.x, psi0g))         stop("'psi.beta' can only be supplied for each group if z.init=list")
  }
  if(beta.x) {
    psi.beta       <- temp.psi <- list(psi.hyper(psi.alpha, cov.mat))
  } else {
    psi.beta       <- len.check(psi.beta, psi0g)
  }
  if(mu0.x)  {
    mu.zero        <- mu
  } else {
    mu.zero        <- len.check(mu.zero, mu0g)
  }
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(ad1.x) {
      alpha.d1     <- list(2)
    } else {
      alpha.d1     <- len.check(alpha.d1, delta0g, P.dim=FALSE)
    }
    if(adk.x) {
      alpha.dk     <- list(10)
    } else {
      alpha.dk     <- len.check(alpha.dk, delta0g, P.dim=FALSE)
    }
  }
  if(!is.element(method, c("FA", "IFA"))) {
    if(verbose)                     cat(paste0("Initialising...\n"))
    clust          <- list()
    pi.alpha       <- list()
    pi.prop        <- list()
    zi             <- list()
    for(g in seq_along(range.G)) {
      G            <- range.G[g]
      if(z.init    == "kmeans")     {
        k.res      <- try(kmeans(dat, G, nstart=100), silent=TRUE)
        if(!inherits(k.res, "try-error"))  {
          zi[[g]]  <- as.numeric(factor(k.res$cluster, levels=seq_len(G)))
        } else                      warning("Cannot initialise cluster labels using kmeans. Try another z.init method", call.=FALSE)
      } else if(z.init  == "list")   {
        zi[[g]]    <- as.numeric(z.list[[g]])
      } else if(z.init  == "mclust") {
        m.res      <- try(Mclust(dat, G), silent=TRUE)
        if(!inherits(m.res, "try_error"))  {
          zi[[g]]  <- as.numeric(m.res$classification)
        } else                      warning("Cannot initialise cluster labels using mclust. Try another z.init method", call.=FALSE)
      } else {
        zips       <- rep(1, N)
        if(!is.element(method, c("IMFA", "IMIFA"))) {
          while(all(length(unique(zips)) != G,
                any(prop.table(tabulate(zips, nbins=G)) < 1/G^2))) {
            pies   <- sim.pi(pi.alpha=rep(alpha.pi, G), nn=0)
            zips   <- sim.z.p(N=N, prob.z=pies)
          }  
        } else {
          if(alpha.pi <= 1)         warning(paste0("Suggestion: increase alpha.pi from ", alpha.pi, " if initialising from the stick-breaking prior"))
          pies     <- sim.pi(pi.alpha=alpha.pi, nn=rep(0, trunc.G), inf.G=TRUE)
          zips     <- sim.z.p(N=N, prob.z=pies)
        }
        zi[[g]]    <- as.numeric(zips)
        rm(zips)
      }
      pi.prop[[g]] <- t(prop.table(tabulate(zi[[g]], nbins=G)))
      mu[[g]]      <- do.call(cbind, lapply(seq_len(G), function(gg) if(pi.prop[[g]][,gg] > 0) colMeans(dat[zi[[g]] == gg,, drop=FALSE]) else rep(0, P)))
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]] else do.call(cbind, lapply(seq_len(G), function(gg) colMeans(dat)))  
      }
      if(beta.x)  {
        if(psi0g) {
          cov.gg   <- lapply(seq_len(G), function(gg) if(pi.prop[[g]][,gg] > 0) cov(dat[zi[[g]] == gg,, drop=FALSE]) else cov.mat)
          psi.beta[[g]] <- do.call(cbind, lapply(seq_len(G), function(gg) psi.hyper(psi.alpha, cov.gg[[gg]])))
        } else {
          psi.beta[[g]] <- replicate(G, temp.psi[[1]])
        }
      }
      if(ad1.x)   {
        alpha.d1[[g]]   <- rep(unlist(alpha.d1), G)
      }
      if(adk.x)   {
        alpha.dk[[g]]   <- rep(unlist(alpha.dk), G)
      }
      clust[[g]]   <- list(z = zi[[g]], pi.alpha = alpha.pi, pi.prop = pi.prop[[g]])
      if(is.element(method, c("MFA", "MIFA"))) {
        sw0g.tmp   <- sw0gs
        if(all(g > 9, any(sw0gs))) {
          sw0g.tmp <- setNames(rep(F, 4), names(sw0gs))
                                    warning(paste0(names(which(sw0gs)), " set to FALSE where G > 9, as 'exact' label-switching is not possible in this case\n"), call.=FALSE)
        }
        clust[[g]] <- append(clust[[g]], list(l.switch = sw0g.tmp))
      }
      if(method == "MIFA") {
        clust[[g]] <- append(clust[[g]], list(alpha.d1 = alpha.d1[[g]], alpha.dk = alpha.dk[[g]]))
      }
    }
  }
  if(is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) {
    mu.zero        <- list(mu.zero[[1]][,1])
    psi.beta       <- list(psi.beta[[1]][,1])
    if(!is.element(method, c("OMFA", "IMFA"))) {
      alpha.d1     <- list(alpha.d1[[1]][1])
      alpha.dk     <- list(alpha.dk[[1]][1])
    }
  } 
  if(all(round(unlist(sapply(mu.zero, sum))) == 0)) {
    mu.zero        <- lapply(mu.zero, function(x) 0)
  }
  if(any(is.na(unlist(psi.beta)))) {
    psi.beta       <- lapply(psi.beta, function(x) replace(x, is.na(x), 0))
  }
  if(any(unlist(psi.beta)   <= 0))  stop("'psi.beta' must be strictly positive")
  if(is.element(method, c("IFA", "MIFA", "IMIFA", "OMIFA"))) {
    if(any(unlist(alpha.d1) <= 0))  stop("'alpha.d1' must be strictly positive")
    if(any((unlist(alpha.dk) 
                   <= beta.dk)))    stop(paste0(ifelse(delta0g, "all 'alpha.dk' values", "alpha.dk"), " must be greater than 'beta.dk'=", beta.dk))
    deltas         <- lapply(seq_along(range.G), function(g) list(alpha.d1 = alpha.d1[[g]], alpha.dk = alpha.dk[[g]]))
  }
  init.time        <- proc.time() - init.start
  
  if(profile)  Rprof()
  if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA"))) {
    if(length(range.G) == 1)  {
      start.time   <- proc.time()
      if(meth[Gi]  != "MIFA") {
        gibbs.arg  <- append(temp.args, deltas[[Gi]])
      }
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),                          
                                     args=append(list(data = dat, N = N, G = range.G, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "IFA") clust[[Gi]]), gibbs.arg))
    } else {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        if(meth[Gi]  == "IFA") {
          gibbs.arg  <- append(temp.args, deltas[[Gi]])
        }
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MIFA") clust[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Gi/length(range.G) * 100, 2), "% Complete\n"))
      }
    }
  } else if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")))   {
    if(all(length(range.G) == 1, length(range.Q) == 1)) {
      start.time   <- proc.time()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]), 
                                     args=append(list(data = dat, N = N, G = range.G, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
    } else if(length(range.G) == 1) {
      start.time   <- proc.time()
      for(q in range.Q) { 
        Qi         <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = range.G, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Qi/length(range.Q) * 100, 2), "% Complete\n"))
      }
    } else if(length(range.Q) == 1) {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
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
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
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
    if(length(Labels) != N)         stop(paste0("Labels must be a factor of length N=",  N))
    range.G        <- nlevels(Labels)
    start.time     <- proc.time()
    for(g in seq_len(range.G)) {
      temp.dat     <- dat[Labels == levels(Labels)[g],]
      imifa[[g]]          <- list()
      imifa[[g]][[Qi]]    <- do.call(paste0("gibbs.", "IFA"),
                                     args=append(list(data = temp.dat, N = nrow(temp.dat), Q = range.Q), gibbs.arg))
      if(verbose)                   cat(paste0(round(g/range.G * 100, 2), "% Complete\n"))
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "MFA", length(range.G) * length(range.Q),
                            ifelse(method == "MIFA", length(range.G),
                              ifelse(method == "classify", nlevels(Labels), 
                                     length(range.Q))))
  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }
  dat.name  <- as.character(match.call()$dat)
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")))   {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  } else {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], "IFA"))
  }
  if(!is.element(method, c("FA", "IFA"))) {
    for(g in seq_along(range.G)) {
      attr(imifa[[g]], 
           "Z.init")      <- factor(zi[[g]], levels=seq_len(range.G[g]))
    }
  }
  gnames    <- paste0(range.G, ifelse(range.G == 1, "Group", "Groups"))
  names(imifa)            <- gnames
  attr(imifa, "Center")   <- centering
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors")  <- range.Q
  attr(imifa, "Groups")   <- range.G
  if(!is.element(method, c("FA", "IFA"))) {
    attr(imifa, "Init.Z") <- z.init
    attr(imifa, 
         "Label.Switch")  <- any(sw0gs)
  }
  method                  <- names(table(meth)[max(table(meth))])
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.name
  attr(imifa, "Obs")      <- N
  attr(imifa, "Scaling")  <- scal
  attr(attr(imifa,
  "Scaling"), "Method")   <- scaling
  if(is.element(method, c("IMFA", "IMIFA"))) {
    attr(imifa, 
         "MH.step")       <- MH.step
    attr(imifa,
         "Gen.Slice")     <- gen.slice
  }
  attr(imifa, "Store")    <- length(iters)
  switches                <- c(switches, a.sw = ifelse(is.element(method, c("IMIFA", "IMFA")), MH.step, F))
  if(is.element(method, c("FA", "IFA"))) {
    switches["pi.sw"]     <- F
  }
  attr(imifa, "Switch")   <- switches
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    attr(imifa, "Time")   <- lapply(list(Total = tot.time, Average = avg.time, Z.Initialisation = init.time), function(x) round(x, 2)) 
  } else {
    attr(imifa, "Time")   <- lapply(list(Total = tot.time, Average = avg.time), function(x) round(x, 2)) 
  }
  if(all(length(range.G)  == 1,
         length(range.Q)  == 1)) {
    attr(imifa, "Time")   <- attr(imifa, "Time")[-2]
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