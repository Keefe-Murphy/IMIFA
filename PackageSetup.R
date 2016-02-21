#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages    <- c("MCMCpack", "slam")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  invisible(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  invisible(lapply(setdiff(packages, (.packages())), library, ch=T))
}
rm(packages)

imifa       <- function(dat = NULL, method = c("IMIFA", "MIFA", "MFA", "IFA", "FA", "classify"), n.iters = 50000,
                        Label = NULL, factanal = F, Q.star = NULL, range.Q = NULL, Q.fac = NULL, thinning = 2,
                        burnin = n.iters/5, n.store = ceiling((n.iters - burnin)/thinning),
                        centering = F, scaling = c("unit", "pareto", "none"), verbose = F, adapt = T, b0 = NULL, b1 = NULL, 
                        prop = NULL, epsilon = NULL, sigma.mu = NULL, sigma.l = NULL, psi.alpha = NULL, psi.beta = NULL,
                        phi.nu = NULL, alpha.d1 = NULL, alpha.d2 = NULL, profile = F, 
                        mu.switch = T, f.switch = T, load.switch = T, psi.switch = T, ...) {
  
  method    <- match.arg(method)
  scaling   <- match.arg(scaling)
  if(missing(dat))                  stop("Dataset must be supplied")
  if(!exists(deparse(substitute(dat)),
             envir=.GlobalEnv))     stop(paste0("Object ", match.call()$dat, " not found"))
  if(!is.logical(factanal))         stop("factanal  must be TRUE or FALSE")
  if(!is.logical(centering))        stop("centering must be TRUE or FALSE")
  if(!is.logical(verbose))          stop("verbose   must be TRUE or FALSE")
  if(!is.logical(profile))          stop("profile   must be TRUE or FALSE")
  if(is.element(method, 
     c("FA", "IFA"))     &&
     missing(centering)  && !centering) {
    centering  <- T
  } else if(!is.element(method,
     c("FA", "IFA"))) {
     centering <- F
  }
  if(all(centering, missing(mu.switch))) {
    mu.switch  <- F 
  }
  if(!missing(mu.switch) && !mu.switch   && !centering) {
    mu.switch  <- T                 
                                    warning("Means were stored since centering was not applied")
  }
  if(method == "FA") {
    if(missing(range.Q))            stop("Arg. range.Q must be specified")
    if((!load.switch     || !psi.switch) && length(range.Q) > 1) {
      load.switch <- T
      psi.switch  <- T
                                    warning("Loadings and Uniquenesses were stored since the optimum Q from the range supplied needs to be determined")
    }
  }
  switches  <- c(mu.sw=mu.switch, f.sw=f.switch, l.sw=load.switch, p.sw=psi.switch)
  if(!is.logical(switches))         stop("All logical switches must be TRUE or FALSE")
  
  # Remove non-numeric columns & apply centering & scaling if necessary 
  dat       <- as.data.frame(dat)
  dat       <- dat[sapply(dat, is.numeric)]
  if(scaling == "pareto") {
    scaling <- sqrt(as.matrix(apply(dat, 2, sd)))
  } else if(scaling == "unit") {
    scaling <- T
  } else {
    scaling <- F
  }
  dat       <- scale(dat, center=centering, scale=scaling)
  
  # Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  N         <- nrow(dat)
  P         <- ncol(dat)
  if(is.null(rownames(dat))) rownames(dat) <- c(1:N)
  if(missing("sigma.mu"))    sigma.mu      <- 0.5
  if(missing("psi.alpha"))   psi.alpha     <- 4
  if(missing("psi.beta"))    psi.beta      <- 1
  if(method == "FA") {
    if(missing("sigma.l"))   sigma.l       <- 0.5
  } else if(method == "IFA" ||
            method == "classify") {
    if(missing("phi.nu"))    phi.nu        <- 3
    if(missing("alpha.d1"))  alpha.d1      <- 2
    if(missing("alpha.d2"))  alpha.d2      <- 10
    if(missing("b0"))        b0            <- 0.1
    if(missing("b1"))        b1            <- 0.00005
    if(missing("prop"))      prop          <- 3/4
    if(missing("epsilon"))   epsilon       <- ifelse(centering, 0.1, 0.01)
  }
  if(method == "classify") {
    source(paste(getwd(), "/IMIFA-GIT/FullConditionals_", "IFA", ".R", sep=""), local=T)
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", "IFA", ".R", sep=""), local=T)
  } else {
    source(paste(getwd(), "/IMIFA-GIT/FullConditionals_", method, ".R", sep=""), local=T)
    source(paste(getwd(), "/IMIFA-GIT/Gibbs_", method, ".R", sep=""), local=T)
  }
  gibbs.arg <- list(n.iters=n.iters, P=P, sigma.mu=sigma.mu, 
                    psi.alpha=psi.alpha, psi.beta=psi.beta, burnin=burnin, 
                    thinning=thinning, n.store=n.store, verbose=verbose, sw=switches)
  if(profile)  Rprof()
  if(method == "IFA" ||
     method == "classify") {
     gibbs.arg     <- append(gibbs.arg, list(phi.nu=phi.nu, alpha.d1=alpha.d1, alpha.d2=alpha.d2,
                                             adapt=adapt, b0=b0, b1=b1, prop=prop, epsilon=epsilon))
  if(missing(Q.star)) {
     Q.star        <- min(floor(3 * log(P)), P)
  } else if(Q.star  > P)            stop("Number of factors must be less than the number of variables")
    if(!is.logical(adapt))          stop("Arg. must be TRUE or FALSE")
    if(method == "IFA") {
      imifa        <- vector("list", length(Q.star))
      start.time   <- proc.time()
      imifa[[1]]   <- do.call(paste0("gibbs.", method),                          
                              args=append(list(data=dat, N=N, Q=Q.star), gibbs.arg))
    } else {
      if(missing(Label))            stop("Data must be labelled for classification")
      if(!exists(deparse(substitute(Label)),
                 envir=.GlobalEnv)) stop(paste0("Object ", match.call()$Label, " not found"))
      Label   <- as.factor(Label)
      if(length(Label) != N)        stop(paste0("Labels must be a factor of length N=",  n.obs))
      imifa   <- vector("list", nlevels(Label))
      start.time   <- proc.time()
      for(i in 1:nlevels(Label)) {
        temp.dat   <- dat[Label == levels(Label)[i],]
        imifa[[i]] <- do.call(paste0("gibbs.", "IFA"),
                              args=append(list(data=temp.dat, N=nrow(temp.dat), Q=Q.star), gibbs.arg))
        if(verbose)                 cat(paste0(round(i/nlevels(Label) * 100, 2), "% Complete\n"))
      }
    }
  } else if(method == "FA") {
    if((length(range.Q)  == 1 && range.Q >= P) || 
       (length(range.Q)   > 1 && any(range.Q) >= P))   
                                    stop ("Number of factors must be less than the number of variables")
    gibbs.arg      <- append(gibbs.arg, list(sigma.l))
    imifa          <- vector("list", length(range.Q))
    if(length(range.Q)   == 1) {
      start.time   <- proc.time()
      imifa[[1]]   <- do.call(paste0("gibbs.", method), 
                              args=append(list(data=dat, N=N, Q=range.Q), gibbs.arg))
    } else {
      start.time   <- proc.time()
      for(q in range.Q) { 
        Q.ind      <- q - min(range.Q) + 1
        imifa[[Q.ind]]   <- do.call(paste0("gibbs.", method),
                                    args=append(list(data=dat, N=N, Q=q), gibbs.arg))
        if(verbose)                 cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete\n"))
      }
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "FA", length(range.Q), ifelse(method == "classify", nlevels(Label), length(Q.star)))
  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }
    
  dat.name  <- as.character(match.call()$dat)
  attr(imifa, "Center")  <- centering
  attr(imifa, "Date")    <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors") <- if(method == "FA") range.Q else Q.star
  attr(imifa, "Method")  <- paste0(toupper(substr(method, 1, 1)),
                                   substr(method, 2, nchar(method)))
  attr(imifa, "Name")    <- dat.name
  attr(imifa, "Obs")     <- N
  attr(imifa, "Scaling") <- scaling
  attr(imifa, "Store")   <- n.store
  attr(imifa, "Switch")  <- switches
  if(method == "IFA" ||
     length(range.Q) == 1) {
    attr(imifa, "Time")  <- tot.time
  } else {
    attr(imifa, "Time")  <- list(Total = tot.time, Average = avg.time) 
  }
  attr(imifa, "Vars")    <- P
  if(verbose)               print(attr(imifa, "Time"))  
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac)) factanal <- T
  if(factanal) {
    if(missing(Q.fac)) Q.fac   <- round(sqrt(P))
    fac     <- factanal(dat, factors=Q.fac, control=list(nstart=50))
    imifa   <- append(imifa, list(fac = fac))
    names(imifa)         <- c(paste0("IMIFA", 1:(length(imifa) - 1)), "fac")
  } 
  class(imifa)           <- "IMIFA"
  return(imifa)
}

source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/Simulate_Data.R", sep=""))