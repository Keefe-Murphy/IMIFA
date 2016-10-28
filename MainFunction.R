########################################################
### Main Function for Keefe Murphy's IMIFA R Package ###
########################################################

mcmc.IMIFA  <- function(dat = NULL, method = c("IMIFA", "IMFA", "OMIFA", "OMFA", "MIFA", "MFA", "IFA", "FA", "classify"), n.iters = 25000, range.G = NULL, range.Q = NULL,
                        burnin = n.iters/5, thinning = 2, centering = TRUE, scaling = c("unit", "pareto", "none"), mu.zero = NULL, sigma.mu = NULL, sigma.l = NULL, alpha = NULL, 
                        z.list = NULL, z.init = c("kmeans", "list", "mclust", "priors"), psi.alpha = NULL, psi.beta = NULL, adapt = TRUE, b0 = NULL, b1 = NULL, prop = NULL, epsilon = NULL, 
                        nu = NULL,  alpha.d1 = NULL, alpha.d2 = NULL, adapt.at = NULL, beta.d1 = NULL, beta.d2 = NULL, alpha.step = c("gibbs", "metropolis", "fixed"), alpha.hyper = NULL, 
                        ind.slice = TRUE, rho = NULL, DP.lab.sw = TRUE, trunc.G = NULL, profile = FALSE, verbose = TRUE, discount = NULL, learn.d = FALSE, d.hyper = NULL, mu0g = FALSE, 
                        psi0g = FALSE, delta0g = FALSE, mu.switch = TRUE, score.switch = TRUE, load.switch = TRUE, psi.switch = TRUE, pi.switch = TRUE, factanal = FALSE, Q.fac = NULL) {
  
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
  dat.nam   <- gsub("[[:space:]]", "", deparse(substitute(dat)))
  nam.dat   <- gsub("\\[.*", "", dat.nam)
  pattern   <- c("(", ")")
  nam.x     <- gsub(".*\\[(.*)\\].*", "\\1)", dat.nam)
  if(!exists(nam.dat,
     envir=.GlobalEnv))             stop(paste0("Object ", match.call()$dat, " not found"))
  if(any(unlist(vapply(seq_along(pattern), function(p) grepl(pattern[p], nam.dat, fixed=TRUE), logical(1))), 
         !identical(dat.nam, nam.dat) && (any(grepl("[[:alpha:]]", gsub('c', '',  nam.x))) || grepl(":", 
         nam.x,    fixed=TRUE))))   stop("Extremely inadvisable to supply 'dat' subsetted by any means other than row/column numbers or c() indexing: best to create new data object")
  zin.miss  <- missing(z.init)
  zli.miss  <- missing(z.list)
  if(!zli.miss) {
    z.nam   <- gsub("[[:space:]]", "", deparse(substitute(z.list)))
    nam.z   <- gsub("\\[.*", "", z.nam)
    nam.zx  <- gsub(".*\\[(.*)\\].*", "\\1)",   z.nam)
    if(!exists(nam.z,
               envir=.GlobalEnv))   stop(paste0("Object ", match.call()$z.list, " not found"))
    if(any(unlist(vapply(seq_along(pattern), function(p) grepl(pattern[p], nam.z,   fixed=TRUE), logical(1))), 
           !identical(z.nam,   nam.z) && (any(grepl("[[:alpha:]]", gsub('c', '', nam.zx))) || grepl(":",
           nam.zx, fixed=TRUE))))   stop("Extremely inadvisable to supply 'z.list' subsetted by any means other than row/column numbers or c() indexing: best to create new object")
    if(!is.list(z.list))     z.list        <- lapply(list(z.list), as.factor)
    if(!zin.miss &&
        z.init  != "list") { z.init        <- "list"
                                    message("'z.init' set to 'list' as 'z.list' was supplied") 
    }                
  }
  if(!is.logical(factanal))         stop("'factanal' must be TRUE or FALSE")
  if(!is.logical(centering))        stop("'centering' must be TRUE or FALSE")
  if(!is.logical(verbose))          stop("'verbose' must be TRUE or FALSE")
  if(!is.logical(profile))          stop("'profile' must be TRUE or FALSE")
  
# Remove non-numeric columns & apply centering & scaling if necessary 
  burnin    <- as.integer(burnin)
  thinning  <- as.integer(thinning)
  if(!is.numeric(n.iters))          stop("'n.iters' must be numeric")
  n.iters   <- max(burnin + 1, as.integer(n.iters))
  iters     <- seq(from=burnin + 1, to=n.iters, by=thinning)
  iters     <- iters[iters > 0]
  raw.dat   <- as.data.frame(dat)
  num.check <- vapply(raw.dat, is.numeric, logical(1))
  if(sum(num.check) != ncol(dat)) { message("Non-numeric columns removed")
    raw.dat <- raw.dat[num.check]
  }
  if(length(iters)  <= 1)           stop("Run a longer chain!")
  if(any(is.na(raw.dat)))  {        message("Rows with missing values removed from data")
    raw.dat <- raw.dat[complete.cases(raw.dat),]
  }          
  if(method != "classify") {
    if(scaling   != "none")   {
      scal  <- apply(raw.dat, 2, sd)
      if(scaling == "pareto") {
       scal <- sqrt(scal)
      }
    } else {
      scal  <- FALSE
    }
    dat     <- scale(raw.dat, center=centering, scale=scal)
  } else   {
    dat     <- raw.dat
  }
  centered  <- if(method == "classify") round(colSums(dat)) == 0 else any(centering, round(colSums(dat)) == 0)
  N         <- nrow(dat)
  P         <- ncol(dat)
  lnN       <- log(N)
  
# Manage storage switches & warnings for other function inputs
  if(!missing(mu.switch) && all(!mu.switch, ifelse(method == "classify", 
     !centering, !centered)))       warning("Centering hasn't been applied - are you sure you want mu.switch=FALSE?", call.=FALSE)
  switches  <- c(mu.sw=mu.switch, s.sw=score.switch, l.sw=load.switch, psi.sw=psi.switch, pi.sw=pi.switch)
  if(!is.logical(switches))         stop("All logical switches must be TRUE or FALSE")
  if(N < 2)                         stop("Must have more than one observation")
  G.x       <- missing(range.G)
  alpha.x   <- missing(alpha.step)
  alpha.step       <- match.arg(alpha.step)
  if(!is.logical(learn.d))          stop("'learn.d' must be TRUE or FALSE")
  if(learn.d)                       stop("Pitman-Yor discount hyperparameter must remain fixed; learning not yet implemented")
  if(missing(d.hyper))       d.hyper       <- c(1, 1)
  if(any(d.hyper   <= 0))           stop("'Discount Beta prior hyperparameters must be strictly positive")
  discount         <- ifelse(missing(discount), ifelse(learn.d, rbeta(1, d.hyper[1], d.hyper[2]), 0), discount)
  if(discount       < 0 || 
     discount      >= 1)            stop("'discount' must lie in the interval [0, 1)")
  if(all(!is.element(method, c("IMFA", "IMIFA")), alpha.step != "fixed"))  {
    alpha.step     <- "fixed"
    if(!alpha.x)                    warning(paste0("'alpha.step' must be given as 'fixed' for the ", method, " method"), call.=FALSE)
  }
  if(!is.element(method, c("MFA", "MIFA")))      {
    if(length(range.G) > 1)         stop(paste0("Only one 'range.G' value can be specified for the ", method, " method"))
    if(all(!G.x, is.element(method, c("FA", "IFA"))) &&  
       range.G  > 1)                warning(paste0("'range.G' must be 1 for the ", method, " method"), call.=FALSE)
    if(is.element(method, c("OMIFA", "OMFA", "IMFA", "IMIFA"))) {
      if(G.x) {
        range.G    <- ifelse(N <= 51, N - 1, max(20, ceiling(3 * log(N))))
      }
      if(range.G    < ceiling(lnN)) stop(paste0("'range.G' should be at least log(N) (=log(", N, "))", " for the ", method, " method"))
      if(is.element(method, c("IMFA", "IMIFA"))) {
        if(!is.logical(ind.slice))  stop("'ind.slice' must be TRUE or FALSE")
        if(!is.logical(DP.lab.sw))  stop("'DP.lab.sw' must be TRUE or FALSE")
        if(missing(rho)) {
          rho      <- 0.75
        }
        if(all(length(rho) > 1,
           rho > 1 && rho <= 0))    stop("'rho' must be a single number in the interval (0, 1]")
        if(missing(alpha.hyper))    {
          alpha.hyper     <- if(alpha.step == "gibbs") c(2, 1) else if(alpha.step == "metropolis") c(- discount, range.G/2) else c(0, 0)
        }
        if(discount > 0) {
          alpha.hyper     <- shift.gamma(shape=alpha.hyper[1], rate=alpha.hyper[2], shift=discount)
        }
        a.len      <- length(alpha.hyper)
        if(a.len   != 2)            stop(paste0("'alpha.hyper' must be a vector of length 2, giving the ", ifelse(alpha.step == "gibbs", "shape and rate hyperparameters of the gamma prior for alpha when alpha.step is given as 'gibbs'", ifelse(alpha.step == "metropolis", "lower and upper limits of the uniform prior/proposal for alpha when alpha.step is given as 'metropolis'")))) 
        a.hyp.1    <- alpha.hyper[1]
        a.hyp.2    <- alpha.hyper[2]
        if(alpha.step == "gibbs")   {
          if(a.hyp.1  <= 0)         stop("The shape of the gamma prior for alpha must be strictly positive")
          if(a.hyp.2  <= 0)         stop("The rate of the gamma prior for alpha must be strictly positive")
        }
        if(alpha.step == "metropolis") {
          if(a.hyp.1   < -discount) stop(paste0("The lower limit of the uniform prior/proposal for alpha must be ", ifelse(discount == 0, "strictly positive", paste0("greater than -discount (=", - discount, ")"))))
          if(a.hyp.2   < 1)         stop("The upper limit of the uniform prior/proposal for alpha must be at least 1")
          if(a.hyp.2  <= a.hyp.1)   stop(paste0("The upper limit (=", a.hyp.2, ") of the uniform prior/proposal for alpha must be greater than the lower limit (=", a.hyp.1, ")"))  
        }
        if(missing(trunc.G))  {
          trunc.G  <- ifelse(N < 100, N, 100)
        }
        if(length(trunc.G) > 1)     stop("'trunc.G' must be a single number")
        if(ifelse(N > 100, trunc.G < 100,
                  trunc.G  < N))    warning(paste0("'trunc.G' should only be less than min(N=", N, ", 100) for practical reasons in heavy computational/memory burden cases"), call.=FALSE)
        if(trunc.G  < range.G)      stop(paste0("'trunc.G' must be at least range.G=", range.G))
        if(trunc.G  > N)            stop(paste0("'trunc.G' cannot be greater than N=", N))
      }
    } else if(method == "classify") {
      if(!zin.miss && 
         z.init    != "list") {     stop("'z.init' must be set to 'list' for classification")
      } else z.init       <- "list"
      if(zli.miss)                  stop("Data labels must be supplied via 'z.list' for classification")
      levs         <- nlevels(z.list)
      if(length(z.list)    > 1)     stop("Only one set of labels can be supplied via 'z.list'")
      zlabels      <- unlist(z.list)
      if(length(zlabels)  != N)     stop(paste0("'z.list' must be a factor of length N=",  N)) 
      if(!missing(range.G) && any(length(range.G > 1), 
          range.G  != levs))   {    warning("Forced 'range.G' equal to the number of levels in 'zlabels' for the 'classify' method")
      }
      range.G      <- levs
    } else {
      range.G      <- 1
    }
    meth    <- method
  } else {
    alp3    <- 3 * alpha
    if(G.x)                         stop("'range.G' must be specified")
    if(any(range.G  < 1))           stop("'range.G' must be strictly positive")
    if(any(range.G  > alp3 * lnN))  warning(paste0("'range.G' MUCH greater than log(N) (=log(", N, ")):\n Empty clusters are likely, consider running an overfitted or infinite mixture"), call.=FALSE)
    range.G <- sort(unique(range.G))
    meth    <- rep(method, length(range.G))                               
  }
  if(any(range.G >= N))             stop(paste0("'range.G' must be less than the number of observations N=", N))
  if(range.G[1]  == 1)   {
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
  Q.min     <- min(ceiling(log(P)), ceiling(log(N)))
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(Q.miss)                      stop("'range.Q' must be specified") 
    if(any(range.Q < 0))            stop(paste0("'range.Q' must be non-negative for the ", method, " method"))
  } else {
    if(Q.miss)        range.Q    <- min(ifelse(P > 500, 12 + floor(log(P)), floor(3 * log(P))), N - 1)
    if(!is.logical(adapt))          stop("'adapt' must be TRUE or FALSE") 
    if(length(range.Q) > 1)         stop(paste0("Only one starting value for 'range.Q' can be supplied for the ", method, " method"))
    if(range.Q    <= 0)             stop(paste0("'range.Q' must be strictly positive for the ", method, " method"))
    if(all(adapt, range.Q < Q.min)) stop(paste0("'range.Q' must be at least min(log(P), log(N)) for the ", method, " method when 'adapt' is TRUE"))
  }
  range.Q   <- sort(unique(range.Q)) 
  len.G     <- ifelse(method == "classify", range.G, length(range.G))
  len.Q     <- length(range.Q)
  len.X     <- len.G * len.Q
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    if(missing("sigma.l"))   sigma.l       <- 1
    if(sigma.l <= 0)                stop("'sigma.l' must be strictly positive")            
  } else {            
    if(missing("nu"))        nu            <- 1.5
    if(nu <= 1)                     stop("'nu' must be strictly greater than 1")
    if(missing("beta.d1"))   beta.d1       <- 1
    if(missing("beta.d2"))   beta.d2       <- 1
    if(missing("b0"))        b0            <- 0.1
    if(b0  < 0)                     stop("'b0' must be positive to ensure valid adaptation probability")
    if(missing("b1"))        b1            <- 0.00005
    if(b1 <= 0)                     stop("'b1' must be strictly positive to ensure adaptation probability decreases")
    if(missing("prop"))      prop          <- 3/4
    if(abs(prop - (1 - prop)) < 0)  stop("'prop' must be a single number between 0 and 1")
    if(missing("adapt.at"))  adapt.at      <- ifelse(is.element(method, c("IFA", "MIFA")), burnin, 0)
    if(missing("epsilon"))   epsilon       <- ifelse(any(centered, centering), 0.1, 0.05)
    if(adapt.at < 0 ||
       adapt.at > burnin)           stop("'adapt.at' must be a single number in the interval [0, burnin]")
    if(epsilon <= 0 ||
       epsilon >= 1)                stop("'epsilon' must be a single number in the interval (0, 1)")
  } 
  if(any(range.Q  >= P)) {          
    if(all(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")),
       isTRUE(adapt)))   {          warning(paste0("Starting value for number of factors is not less than the number of variables, ", P), call.=FALSE)
    } else if(any(is.element(method, c("FA", "MFA", "OMFA", "IMFA")),
              all(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA")), 
                  isTRUE(!adapt)))) stop(paste0("Number of factors must be less than the number of variables, ", P))
  } 
  if(any(range.Q  >= N))            stop(paste0("Number of factors must be less than the number of observations, ", N))
  if(any(all(method == "MFA",  any(range.G > 1)) && any(range.Q > 0),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA",
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["l.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Loadings & Uniquenesses not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Loadings not stored: may be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Uniquenesses not stored: will be unable to estimate covariance matrix and compute error metrics", call.=FALSE)
  }
  if(any(all(method == "MFA",  any(range.G > 1)),
         all(method == "MIFA", any(range.G > 1)), is.element(method, c("IMIFA", 
     "IMFA", "OMIFA", "OMFA"))))  {
    if(all(!switches["mu.sw"], 
           !switches["psi.sw"]))  {
                                    warning("Means & Uniquenesses not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["mu.sw"]) { warning("Means not stored: posterior mean estimates won't be available", call.=FALSE)
    } else if(!switches["psi.sw"])  warning("Uniquenesses not stored: posterior mean estimates won't be available", call.=FALSE)
  }
  if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")) && all(range.Q == 0)) {   
    if(all(switches[c("s.sw", "l.sw")]))  {
                                    warning("Scores & Loadings not stored as model has zero factors", call.=FALSE)
    } else if(switches["s.sw"])   { warning("Scores not stored as model has zero factors", call.=FALSE)
    } else if(switches["l.sw"])   { warning("Loadings not stored as model has zero factors", call.=FALSE)
    }                               
    switches[c("s.sw", "l.sw")]  <- FALSE                              
  } else {
    if(all(!switches[c("s.sw", "l.sw")])) { 
                                    warning("Posterior Scores & Loadings won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["s.sw"])  { warning("Posterior Scores won't be available as they're not being stored", call.=FALSE)
    } else if(!switches["l.sw"])  { warning("Posterior Loadings won't be available as they're not being stored", call.=FALSE)
    }
  }
  if(all(is.element(method, c("FA", "IFA")), 
         !zin.miss || 
         !zli.miss))                message(paste0("z does not need to be initialised for the ", method, " method"))
  if(is.element(method, c("MFA", "MIFA", "classify"))) {
    if(!is.logical(mu0g))           stop("'mu0g' must be TRUE or FALSE")
    if(!is.logical(psi0g))          stop("'psi0g' must be TRUE or FALSE")
    if(!is.logical(delta0g))        stop("'delta0g' must be TRUE or FALSE")  
    if(all(method == "MFA",
           delta0g))                stop("'delta0g' cannot be TRUE for the 'MFA' method")
    if(method == "classify") mu0g          <- TRUE
  }
  sw0gs     <- c(mu0g = mu0g, psi0g = psi0g, delta0g = delta0g)
  if(all(!is.element(method, c("MFA", "MIFA", "classify")), 
         any(sw0gs)))               stop(paste0(names(which(sw0gs)), " should be FALSE for the ", method, " method\n"))
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    if(missing("alpha"))     alpha         <- ifelse(is.element(method, c("OMIFA", "OMFA")), 0.5/range.G, 
                                              ifelse(alpha.step == "gibbs", rgamma(1, a.hyp.1, a.hyp.2) - discount,
                                              ifelse(alpha.step == "metropolis", runif(1, a.hyp.1, a.hyp.2), 1)))
    if(length(alpha) != 1)          stop("'alpha' must be specified as a scalar to ensure an exchangeable prior")
    if(alpha <= -discount)          stop(paste0("'alpha' must be ", ifelse(discount != 0, paste0("greater than -discount (i.e. > ", - discount, ")"), "strictly positive")))
    if(all(is.element(method,  c("IMIFA", "IMFA")),
       alpha.step == "fixed"))      warning(paste0("'alpha' fixed at ", alpha, " as it's not being learned via Gibbs/Metropolis-Hastings updates"), call.=FALSE)
    if(all(!is.element(method, c("IMFA", "IMIFA")),
           alpha   > 1))            warning("Are you sure alpha should be greater than 1?", call.=FALSE)
                             z.init        <- match.arg(z.init)
    if(all(is.element(method,  c("OMIFA", "OMFA")), !is.element(z.init, 
       c("list", "kmeans"))))       stop(paste0("'z.init' must be set to 'list' or 'kmeans' for the ", method, " method to ensure all groups are populated at the initialisation stage"))
    if(!zli.miss) {
      if(length(z.list)   != len.G) {
                                    stop(paste0("'z.list' must be a list of length ", len.G))  }
                             list.levels   <- lapply(z.list, nlevels)
      if(!all(list.levels == range.G))            {
        if(!is.element(method, c("IMIFA", 
                       "IMFA", "OMIFA", "OMFA"))) {
                                    stop(paste0("Each element of 'z.list' must have the same number of levels as 'range.G'")) 
        } else                      stop(paste0("Only ", list.levels, " groups are populated according to z.list, but 'range.G' has been set to ", range.G, ":\n  Reset range.G to this value to avoid redunandtly carrying around empty groups or supply a list with ", range.G, " levels"))
      }
      if(!all(lapply(z.list, length) == N)) {
                                    stop(paste0("Each element of 'z.list' must be a vector of length N=", N)) }
    }
    if(all(zli.miss, z.init == "list"))     {
                                    stop(paste0("'z.list' must be supplied if 'z.init' is set to 'list'")) }
  }
  imifa     <- list(list())
  Gi        <- 1
  Qi        <- 1
  gibbs.arg <- list(P = P, sigma.mu = sigma.mu, psi.alpha = psi.alpha, burnin = burnin, 
                    thinning = thinning, iters = iters, verbose = verbose, sw = switches)
  if(is.element(method, c("IMIFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(trunc.G = trunc.G, rho = rho, ind.slice = ind.slice, alpha.step = alpha.step, learn.d = learn.d, 
                                             DP.lab.sw = DP.lab.sw, a.hyper = alpha.hyper, discount = discount, d.hyper = d.hyper))
  }
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    gibbs.arg      <- append(gibbs.arg, list(nu = nu, beta.d1 = beta.d1, beta.d2 = beta.d2, adapt.at = adapt.at,
                                             adapt = adapt, b0 = b0, b1 = b1, prop = prop, epsilon = epsilon))
    temp.args      <- gibbs.arg
  } else {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l = sigma.l))
  }

  init.start       <- proc.time()  
  mu               <- list(colMeans(dat))
  beta.x           <- missing("psi.beta")
  mu0.x            <- missing("mu.zero")
  ad1.x            <- all(missing("alpha.d1"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA", "classify")))
  adk.x            <- all(missing("alpha.d2"), is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA", "classify")))
  if(all(z.init != "list", any(sw0gs))) {
    if(delta0g)                     stop("'delta0g' can only be TRUE if z.init=list\n")
    if(all(!mu0.x, mu0g))           stop("'mu.zero' can only be supplied for each group if z.init=list")
    if(all(!beta.x, psi0g))         stop("'psi.beta' can only be supplied for each group if z.init=list")
  }
  if(beta.x) {
    psi.beta       <- temp.psi <- list(psi.hyper(psi.alpha, cov.mat))
  } else {
    psi.beta       <- len.check(psi.beta, psi0g, method, P, range.G)
  }
  mu.zero          <- if(mu0.x) mu else len.check(mu.zero, mu0g, method, P, range.G)
  if(!is.element(method, c("FA", "MFA", "OMFA", "IMFA"))) {
    alpha.d1       <- if(ad1.x) list(3) else len.check(alpha.d1, delta0g, method, P, range.G, P.dim=FALSE)
    alpha.d2       <- if(adk.x) list(6) else len.check(alpha.d2, delta0g, method, P, range.G, P.dim=FALSE)
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
        k.res      <- try(suppressWarnings(kmeans(dat, centers=G, iter.max=20, nstart=100)), silent=TRUE)
        if(!inherits(k.res, "try-error"))  {
          zi[[g]]  <- as.numeric(factor(k.res$cluster, levels=seq_len(G)))
        } else                      stop("Cannot initialise cluster labels using kmeans. Try another z.init method")
      } else if(z.init  == "list")   {
        zi[[g]]    <- as.numeric(z.list[[g]])
      } else if(z.init  == "mclust") {
        m.res      <- try(Mclust(dat, G), silent=TRUE)
        if(!inherits(m.res, "try_error"))  {
          zi[[g]]  <- as.numeric(m.res$classification)
        } else                      stop("Cannot initialise cluster labels using mclust. Try another z.init method")
      } else {
        zips       <- rep(1, N)
        if(!is.element(method, c("IMFA", "IMIFA"))) {
          while(all(length(unique(zips)) != G,
                any(prop.table(tabulate(zips, nbins=G)) < 1/G^2))) {
            pies   <- sim.pi(pi.alpha=rep(alpha, G), nn=0)
            zips   <- sim.z.p(N=N, prob.z=pies)
          }  
        } else {
          if(alpha <= 1)            warning("Suggestion: supply a value > 1 for 'alpha' if initialising labels from the stick-breaking prior")
          pies     <- sim.pi(pi.alpha=alpha, nn=rep(0, trunc.G), inf.G=TRUE)
          zips     <- sim.z.p(N=N, prob.z=pies)
        }
        zi[[g]]    <- as.numeric(zips)
        rm(zips)
      }
      nngs         <- tabulate(zi[[g]], nbins=G)
      pi.prop[[g]] <- prop.table(nngs)
      mu[[g]]      <- vapply(seq_len(G), function(gg) if(nngs[gg] > 0) colMeans(dat[zi[[g]] == gg,, drop=FALSE]) else rep(0, P), numeric(P))
      if(mu0.x)   {
        mu.zero[[g]]    <- if(mu0g) mu[[g]] else vapply(seq_len(G), function(gg) colMeans(dat), numeric(P))
      }
      if(beta.x)  {
        if(psi0g) {
          cov.gg   <- lapply(seq_len(G), function(gg) if(nngs[gg] > 1) cov(dat[zi[[g]] == gg,, drop=FALSE]) else cov.mat)
          psi.beta[[g]] <- vapply(seq_len(G), function(gg) psi.hyper(psi.alpha, cov.gg[[gg]]), numeric(P))
        } else {
          psi.beta[[g]] <- replicate(G, temp.psi[[1]])
        }
      }
      if(ad1.x)   {
        alpha.d1[[g]]   <- rep(unlist(alpha.d1), G)
      }
      if(adk.x)   {
        alpha.d2[[g]]   <- rep(unlist(alpha.d2), G)
      }
      clust[[g]]   <- list(z = zi[[g]], pi.alpha = alpha, pi.prop = pi.prop[[g]])
      if(is.element(method, c("MFA", "MIFA"))) {
        sw0g.tmp   <- sw0gs
        if(all(g > 9, any(sw0gs))) {
          sw0g.tmp <- setNames(rep(FALSE, 4), names(sw0gs))
                                    warning(paste0(names(which(sw0gs)), " set to FALSE where G > 9, as 'exact' label-switching is not possible in this case\n"), call.=FALSE)
        }
        clust[[g]] <- append(clust[[g]], list(l.switch = sw0g.tmp))
      }
      if(is.element(method, c("classify", "MIFA"))) {
        clust[[g]] <- append(clust[[g]], list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
      }
    }
  }
  if(is.element(method, c("IMIFA", "IMFA", "OMIFA", "OMFA"))) {
    mu.zero        <- list(mu.zero[[1]][,1])
    psi.beta       <- list(psi.beta[[1]][,1])
    if(!is.element(method, c("OMFA", "IMFA"))) {
      alpha.d1     <- list(alpha.d1[[1]][1])
      alpha.d2     <- list(alpha.d2[[1]][1])
    }
  } 
  if(all(round(vapply(mu.zero, sum, numeric(1))) == 0)) {
    mu.zero        <- if(method == "classify") matrix(0, nr=1, nc=range.G) else lapply(mu.zero, function(x) 0)
  }
  if(any(is.na(unlist(psi.beta)))) {
    psi.beta       <- lapply(psi.beta, function(x) replace(x, is.na(x), 0))
  }
  if(any(unlist(psi.beta)   <= 0))  stop("'psi.beta' must be strictly positive")
  if(is.element(method, c("classify", "IFA", "MIFA", "IMIFA", "OMIFA"))) {
    if(!all(MGP.check(unlist(alpha.d1), unlist(alpha.d2), unique(range.Q),
            beta.d1, beta.d2)))     stop("Invalid shrinkage hyperparameter values will not encourage loadings column removal.\n Try using the MGP.check() function in advance.")
    deltas         <- lapply(seq_along(range.G), function(g) list(alpha.d1 = alpha.d1[[g]], alpha.d2 = alpha.d2[[g]]))
  }
  init.time        <- proc.time() - init.start
  fac.time         <- 0
  
  if(profile)  Rprof()
  if(is.element(method, c("IFA", "MIFA", "OMIFA", "IMIFA"))) {
    if(len.G == 1)  {
      start.time   <- proc.time()
      if(meth[Gi]  != "MIFA") {
        gibbs.arg  <- append(temp.args, deltas[[Gi]])
      }
      imifa[[Gi]][[Qi]]   <- do.call(paste0("gibbs.", meth[Gi]),                          
                                     args=append(list(data = dat, N = N, G = range.G, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "IFA") clust[[Gi]]), gibbs.arg))
      fac.time     <- fac.time + imifa[[Gi]][[Qi]]$time
    } else {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        if(meth[Gi]  == "IFA") {
         gibbs.arg <- append(temp.args, deltas[[Gi]])
        }
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] == "MIFA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Gi  != len.G) cat(paste0("Model ", Gi, " of ", len.G, " complete"), "Initialising...", sep="\n")
      }
    }
  } else if(is.element(method, c("FA", "MFA", "OMFA", "IMFA")))   {
    if(all(len.G == 1, len.Q == 1)) {
      start.time   <- proc.time()
      imifa[[Gi]][[Qi]]   <- do.call(paste0("gibbs.", meth[Gi]), 
                                     args=append(list(data = dat, N = N, G = range.G, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
      fac.time     <- fac.time + imifa[[Gi]][[Qi]]$time
    } else if(len.G == 1) {
      start.time   <- proc.time()
      for(q in range.Q)   { 
        Qi         <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = range.G, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Qi  != len.Q) cat(paste0("Model ", Qi, " of ", len.Q, " complete"), "Initialising...", sep="\n")
      }
    } else if(len.Q == 1) {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && Gi  != len.G) cat(paste0("Model ", Gi, " of ", len.G, " complete"), "Initialising...", sep="\n")
      }
    } else {
      mi           <- 0
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        for(q in range.Q) {
          Qi       <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", meth[Gi]),
                                     args=append(list(data = dat, N = N, G = g, Q = q, mu = mu[[Gi]], mu.zero = mu.zero[[Gi]],
                                                      psi.beta = psi.beta[[Gi]], cluster = if(meth[Gi] != "FA") clust[[Gi]]), gibbs.arg))
        mi         <- mi + 1
        fac.time   <- fac.time + imifa[[Gi]][[Qi]]$time
        if(verbose && mi  != len.X) cat(paste0("Model ", mi, " of ", len.X, " complete"), "Initialising...", sep="\n")
        }
      }
    }
  } else if(method == "classify") { stop("'classify' method not yet implemented")   
    start.time     <- proc.time()
    if(centered)                    warning("Data supplied is globally centered, are you sure?", call.=FALSE)
    for(g in seq_len(range.G))  {
      tmp.dat      <- raw.dat[zlabels == levels(zlabels)[g],]
      if(scaling   != "none")   {
        scal  <- apply(tmp.dat, 2, sd)
        if(scaling == "pareto") {
         scal <- sqrt(scal)
        }
      } else {
        scal  <- FALSE
      }
      tmp.dat <- scale(tmp.dat, center=centering, scale=scal)
      imifa[[g]]          <- list()
      gibbs.arg    <- append(temp.args, lapply(deltas[[Gi]], "[[", g))
      imifa[[g]][[Qi]]    <- do.call(paste0("gibbs.", "IFA"),
                                     args=append(list(data = tmp.dat, N = nrow(tmp.dat), mu = mu[[Gi]][,g], mu.zero = mu.zero[[Gi]][,g],
                                                      Q = range.Q, psi.beta = psi.beta[[Gi]][,g]), gibbs.arg))
      fac.time     <- fac.time + imifa[[g]][[Qi]]$time
      if(verbose   && g   != len.G) cat(paste0("Model ", g, " of ", len.G, " complete"), "Initialising...", sep="\n")
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "MFA", len.X,
                           ifelse(method == "MIFA", len.G,
                             ifelse(method == "classify", range.G, len.Q)))
  tot.time  <- tot.time    + init.time
  init.time <- init.time   + fac.time
  for(g in length(imifa)) {
   for(q in length(imifa[[g]])) {
     imifa[[g]][[q]]$time <- NULL
   }
  }

  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }
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
  gnames    <- if(method  == "classify") paste0("Group ", seq_len(range.G)) else paste0(range.G, ifelse(range.G == 1, "Group", "Groups"))
  names(imifa)            <- gnames
  attr(imifa, 
       "Alpha.step")      <- alpha.step
  attr(imifa, "Alpha")    <- if(alpha.step == "fixed") alpha
  if(method == "classify") {
    attr(imifa,
         "Class.Props")   <- tabulate(z.list[[1]], range.G)/N
  }
  attr(imifa, "Center")   <- any(centered, centering)
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa,
       "Disc.step")       <- learn.d
  attr(imifa, "Discount") <- if(!learn.d) discount
  attr(imifa, "Factors")  <- range.Q
  attr(imifa, 
       "Ind.Slice")       <- all(is.element(method, c("IMFA", "IMIFA")), ind.slice)
  attr(imifa, "Groups")   <- range.G
  if(!is.element(method, c("FA", "IFA"))) {
    attr(imifa, "Init.Z") <- z.init
    attr(imifa, 
         "Label.Switch")  <- any(sw0gs)
  }
  method                  <- names(table(meth)[which.max(table(meth))])
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.nam
  attr(imifa, "Obs")      <- N
  attr(imifa, "Scaling")  <- scal
  attr(attr(imifa,
  "Scaling"), "Method")   <- scaling
  attr(imifa, "Store")    <- length(iters)
  switches                <- c(switches, a.sw = alpha.step != "fixed")
  if(is.element(method, c("FA", "IFA"))) {
    switches["pi.sw"]     <- FALSE
  }
  attr(imifa, "Switch")   <- switches
  times                   <- lapply(list(Total = tot.time, Average = avg.time, Initialisation = init.time), function(x) round(x, 2)) 
  if(all(len.G  == 1,
         len.Q  == 1)) {
    times                 <- times[-2]
  }
  if(is.element(method, c("FA", "IFA", "classify"))) {
    times                 <- times[-length(times)]
  }
  class(times)            <- "listof"
  attr(imifa, "Time")     <- times
  attr(imifa, "Vars")     <- P
  if(verbose)                print(attr(imifa, "Time"))  
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac))    factanal <- TRUE
  if(factanal) {
    if(missing(Q.fac)) {
      Q.fac <- if(missing(range.Q)) round(sqrt(P)) else max(1, max(range.Q))
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