#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages    <- c("abind", "MCMCpack", "slam")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  invisible(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  invisible(lapply(setdiff(packages, (.packages())), library, ch=T))
}
rm(packages)

imifa.mcmc  <- function(dat = NULL, method = c("IMIFA", "MIFA", "MFA", "IFA", "FA", "classify"), 
                        n.iters = 50000, Label = NULL, factanal = F, Q.star = NULL, range.G = NULL, 
                        range.Q = NULL, Q.fac = NULL,  burnin = n.iters/5, thinning = 2, centering = F, 
                        scaling = c("unit", "pareto", "none"), verbose = F, adapt = T, b0 = NULL, 
                        b1 = NULL, prop = NULL, epsilon = NULL, sigma.mu = NULL, sigma.l = NULL, 
                        psi.alpha = NULL, psi.beta = NULL, phi.nu = NULL, alpha.d1 = NULL, alpha.d2 = NULL, 
                        alpha.pie = NULL, z.init = c("kmeans", "priors", "list"), z.list = NULL, 
                        profile = F, mu.switch = T, f.switch = T, load.switch = T, psi.switch = T, ...) {
  
  defpar    <- par(no.readonly = T)
  defop     <- options()
  options(warn=1)
  on.exit(suppressWarnings(par(defpar)))
  on.exit(suppressWarnings(options(defop)), add=T)
  if(method == "classification") {
    method  <- "classify"
  }
  method    <- match.arg(method)
  scaling   <- match.arg(scaling)
  if(missing(dat))                  stop("Dataset must be supplied")
  if(!exists(deparse(substitute(dat)),
             envir=.GlobalEnv))     stop(paste0("Object ", match.call()$dat, " not found"))
  if(!is.logical(factanal))         stop("factanal  must be TRUE or FALSE")
  if(!is.logical(centering))        stop("centering must be TRUE or FALSE")
  if(!is.logical(verbose))          stop("verbose   must be TRUE or FALSE")
  if(!is.logical(profile))          stop("profile   must be TRUE or FALSE")
  if(is.element(method, c("FA", "IFA", "classify"))) {
    if(missing(centering) || !centering) {
     centering <- T  
    } 
  } else {
     centering <- F
  }
  
# Remove non-numeric columns & apply centering & scaling if necessary 
  n.store   <- ceiling((n.iters - burnin)/thinning)
  dat       <- dat[sapply(dat, is.numeric)]
  if(scaling == "pareto") {
    scaling <- sqrt(as.matrix(apply(dat, 2, sd)))
  } else {
    scaling <- scaling == "unit"
  }
  dat       <- as.data.frame(scale(dat, center=centering, scale=scaling))
  N         <- nrow(dat)
  P         <- ncol(dat)
  
# Manage storage switches and warnings for other function inputs
  if(all(centering, missing(mu.switch))) {
    mu.switch  <- F 
  }
  if(!missing(mu.switch) && all(!mu.switch, !centering)) {
    mu.switch  <- T                 
                                    warning("Means were stored since centering was not applied", call.=F)
  }
  switches  <- c(mu.sw=mu.switch, f.sw=f.switch, l.sw=load.switch, p.sw=psi.switch)
  if(!is.logical(switches))         stop("All logical switches must be TRUE or FALSE")
  if(!is.element(method, c("MFA", "MIFA"))) {
    if(!missing(range.G) &&  
       any(range.G  > 1))           warning(paste0("range.G must be 1 for method = ", method), call.=F)
    range.G <- 1
  } else {
    if(missing(range.G))            stop("range.G must be specified")
    if(any(range.G  < 1))           stop("range.G must be strictly positive")
    if(all(range.G == 1)) {
      if(method != "IMIFA") {
                                    warning(paste0("Forced use of ", method, " method as range.G is equal to 1"), call.=F)                        
        if(method  == "MFA") {
          method   <- "FA"                              
        } else if(method == "MIFA") {
          method   <- "IFA"
        } 
      }
    }                               
    range.G <- sort(unique(range.G))
  } 
  no.fac    <- all(all(range.Q == 0), is.element(method, c("FA", "MFA")))
  if(is.element(method, c("FA", "MFA"))) {
    if(!missing(Q.star))            rm(Q.star)
    if(missing(range.Q))            stop("range.Q must be specified")
    if(any(range.Q < 0))            stop("range.Q must be strictly non-negative")
    range.Q <- sort(unique(range.Q))   
  } else {
    if(!missing(range.Q))           rm(range.Q)
    if(missing(Q.star)) {
      Q.star       <- min(floor(3 * log(P)), P, N - 1)
    } else {
      if(Q.star     > P)            stop("Number of factors must be less than the number of variables")
      if(Q.star    >= N)            stop("Number of factors must be less than the number of observations")
    } 
    if(!is.logical(adapt))          stop("adapt must be TRUE or FALSE") 
  }
  if(no.fac) {   
    if(all(switches[c("f.sw", "l.sw")]))  {
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
  
# Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  if(is.null(rownames(dat))) rownames(dat) <- 1:N
  if(missing("sigma.mu"))    sigma.mu      <- 0.5
  if(missing("psi.alpha"))   psi.alpha     <- 4
  if(missing("psi.beta"))    psi.beta      <- 1
  if(is.element(method, c("FA", "MFA"))) {
    if(missing("sigma.l"))   sigma.l       <- 0.5
  } else {
    if(missing("phi.nu"))    phi.nu        <- 3
    if(missing("alpha.d1"))  alpha.d1      <- 2
    if(missing("alpha.d2"))  alpha.d2      <- 10
    if(missing("b0"))        b0            <- 0.1
    if(missing("b1"))        b1            <- 0.00005
    if(missing("prop"))      prop          <- 3/4
    if(missing("epsilon"))   epsilon       <- ifelse(centering, 0.1, 0.005)
  } 
  if(!is.element(method, c("FA", "IFA"))) {
    if(missing("alpha.pie")) alpha.pie     <- 0.5
                             z.init        <- match.arg(z.init)
    if(!missing(z.list)) {
                             z.list        <- lapply(z.list, as.factor)
      if(z.init != "list")          stop(paste0("z.init must be set to 'list' if z.list is supplied"))
      if(length(z.list)   != length(range.G)) {
                                    stop(paste0("z.list must be a list of length ", length(range.G))) }
      if(!all(lapply(z.list, nlevels)      == range.G)) {
                                    stop(paste0("Each element of z.list must have the same number of levels as range.G")) }
      if(!all(lapply(z.list, length)       == N)) {
                                    stop(paste0("Each element of z.list must be a vector of length N=", N)) }
    }
    if(all(missing(z.list),  z.init == "list")) {
                                    stop(paste0("z.list must be supplied if z.init is set to 'list"))
    }
  }  
  
  if(method == "classify") {
    source(paste(getwd(), "/IMIFA-GIT/FullConditionals_", "IFA", ".R", sep=""), local=T)
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", "IFA", ".R", sep=""), local=T)
  } else {
    source(paste(getwd(), "/IMIFA-GIT/FullConditionals_", method, ".R", sep=""), local=T)
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", method, ".R", sep=""), local=T)
  }

  imifa     <- list(list())
  Gi        <- 1
  Qi        <- 1
  gibbs.arg <- list(n.iters = n.iters, P = P, sigma.mu = sigma.mu, 
                    psi.alpha = psi.alpha, psi.beta = psi.beta, burnin = burnin, 
                    thinning = thinning, n.store = n.store, verbose = verbose, sw = switches)
  if(!is.element(method, c("FA", "MFA"))) {
    gibbs.arg     <- append(gibbs.arg, list(phi.nu = phi.nu, alpha.d1 = alpha.d1, alpha.d2 = alpha.d2,
                                          adapt = adapt, b0 = b0, b1 = b1, prop = prop, epsilon = epsilon))
  } else {
    gibbs.arg     <- append(gibbs.arg, list(sigma.l = sigma.l))
  }
  if(!is.element(method, c("FA", "IFA", "classify"))) {
    gibbs.arg      <- append(gibbs.arg, list(alpha.pie = alpha.pie, z.init = z.init))
  }
  
  if(profile)  Rprof()
  if(method == "IFA") {
    start.time     <- proc.time()
    imifa[[Gi]][[Qi]]   <- do.call(paste0("gibbs.", method),                          
                                   args=append(list(data = dat, N = N, Q = Q.star), gibbs.arg))
  }
  if(method == "classify") {
    if(missing(Label))              stop("Data must be labelled for classification")
    if(!exists(deparse(substitute(Label)),
               envir=.GlobalEnv))   stop(paste0("Object ", match.call()$Label, " not found"))
    Label   <- as.factor(Label)
    if(length(Label) != N)          stop(paste0("Labels must be a factor of length N=",  n.obs))
    range.G        <- nlevels(Label)
    start.time     <- proc.time()
    for(g in 1:range.G) {
      temp.dat     <- dat[Label == levels(Label)[g],]
      imifa[[g]]          <- list()
      imifa[[g]][[Qi]]    <- do.call(paste0("gibbs.", "IFA"),
                                     args=append(list(data = temp.dat, N = nrow(temp.dat), Q = Q.star), gibbs.arg))
      if(verbose)                   cat(paste0(round(g/range.G * 100, 2), "% Complete\n"))
    }
  } 
  if(is.element(method, c("FA", "MFA"))) {
    if(any(all(length(range.Q)   == 1, any(range.Q >= P, range.Q >= N - 1)), 
           all(length(range.Q)    > 1, any(any(range.Q  >= P), any(range.Q  >= N - 1)))))   
                                    stop("Number of factors must be less than the number of variables and number of observations")
    if(all(length(range.G) == 1, length(range.Q) == 1)) {
      start.time   <- proc.time()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", method), 
                                     args=append(list(data = dat, N = N, Q = range.Q), gibbs.arg))
    } else if(length(range.G) == 1) {
      start.time   <- proc.time()
      for(q in range.Q) { 
        Qi         <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", method),
                                     args=append(list(data = dat, N = N, Q = q), gibbs.arg))
        if(verbose)                 cat(paste0(round(Qi/length(range.Q) * 100, 2), "% Complete\n"))
      }
    } else if(length(range.Q) == 1) {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        imifa[[Gi]][[Qi]] <- do.call(paste0("gibbs.", method),
                                     args=append(list(data = dat, N = N, G = g, Q = range.Q,
                                                 if(z.list == "list") z.list = z.list[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round(Gi/length(range.G) * 100, 2), "% Complete\n"))
      }
    } else {
      start.time   <- proc.time()
      for(g in range.G) {
        Gi         <- which(range.G == g)
        imifa[[Gi]]       <- list()
        for(q in range.Q) {
          Qi       <- which(range.Q == q)
        imifa[[Gi]][[Qi]] <- do.call(paste0("Gibbs.", method),
                                     args=append(list(data = dat, N = N, G = g, Q = q,
                                                 if(z.list == "list") z.list = z.list[[Gi]]), gibbs.arg))
        if(verbose)                 cat(paste0(round((Gi * Qi)/(length(range.G) * length(range.Q)) * 100, 2), "% Complete\n"))
        }
      }
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "MFA", length(range.G) * length(range.Q),
                          ifelse(method == "FA",  length(range.Q), 
                            ifelse(method == "MIFA", length(range.G),
                              ifelse(method == "classify", nlevels(Label), 
                                     length(Q.star)))))
  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }  
  dat.name  <- as.character(match.call()$dat)
  if(is.element(method, c("FA", "MFA"))) {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], paste0(range.Q, ifelse(range.Q == 1, "Factor", "Factors"))))
  } else {
    imifa   <- lapply(seq_along(imifa), function(x) setNames(imifa[[x]], "IMIFA"))
  }
  gnames    <- paste0(range.G, ifelse(range.G == 1, "Group", "Groups"))
  names(imifa)            <- gnames
  attr(imifa, "Center")   <- centering
  attr(imifa, "Date")     <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors")  <- if(is.element(method, c("FA", "MFA"))) range.Q else Q.star
  attr(imifa, "Groups")   <- if(method != "IMIFA") range.G else G.star
  attr(imifa, "Method")   <- paste0(toupper(substr(method, 1, 1)),
                                    substr(method, 2, nchar(method)))
  attr(imifa, "Name")     <- dat.name
  attr(imifa, "Obs")      <- N
  attr(imifa, "Scaling")  <- scaling
  attr(imifa, "Store")    <- n.store
  attr(imifa, "Switch")   <- switches
  if(is.element(method, c("IFA", "IMIFA")) || length(range.Q) == 1) {
    attr(imifa, "Time")   <- tot.time
  } else {
    attr(imifa, "Time")   <- list(Total = tot.time, Average = avg.time) 
  }
  attr(imifa, "Vars")     <- P
  if(verbose)                print(attr(imifa, "Time"))  
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac))   factanal <- T
  if(factanal) {
    if(missing(Q.fac)) {
      if(missing(range.Q)) {
        Q.fac      <- round(sqrt(P))
      } else {
        Q.fac      <- max(1, max(range.Q))
      }
    }
    fac     <- factanal(dat, factors=Q.fac, control=list(nstart=50))
    imifa   <- append(imifa, list(fac = fac))
    names(imifa)[length(imifa)] <- "Factanal"
  } 
  class(imifa)            <- "IMIFA"
  return(imifa)
}

source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/SimulateData.R", sep=""))