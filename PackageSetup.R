#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

preamble    <- function(seed=21092015, rem.lib=F, rem.all=F, ...) {
  
  if(!is.element(rem.lib, c(T, F)))   stop("Arg. must be TRUE or FALSE")
  if(!is.element(rem.all, c(T, F)))   stop("Arg. must be TRUE or FALSE")
  set.seed(seed)
  assign("def.par", par(), envir=.GlobalEnv)
  packages  <- c("pgmm", "car", "MCMCpack", "compiler")
  if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
    invisible(install.packages(setdiff(packages, rownames(installed.packages()))))
  }
  if(length(setdiff(packages, (.packages()))) > 0) {
    invisible(lapply(setdiff(packages, (.packages())), library, ch=T))
  }
  
# WARNING: Remove loaded libraries
  if(rem.lib) {
    pkgs    <- names(sessionInfo()$otherPkgs)
    pkgs    <- paste('package:', pkgs, sep = "")
    invisible(lapply(pkgs, detach, ch = T, unload = T, force= T))
  }

# WARNING: Remove everything
  if(rem.all) rm(list = ls(all = TRUE))
}

preamble()

imifa.gibbs <- function(dat=NULL, n.iters=50000, method=c("IMIFA", "MIFA", "IFA", "FA"), 
                        factanal=F, Q.star=NULL, range.Q=NULL, Q.fac=NULL, thinning=2,
                        burnin=n.iters/5 - 1, n.store=ceiling((n.iters - burnin)/thinning),
                        centering=T, scaling=T, print=T, adapt=T,
                        sigma.mu=NULL, sigma.l=NULL, psi.alpha=NULL, psi.beta=NULL,
                        phi.nu=NULL, delta.a1=NULL, delta.a2=NULL, profile=F, ...) {
  
  method    <- match.arg(method)
  assign("method", method, envir=.GlobalEnv)
  if(missing(dat))                    stop("Dataset must be supplied")
  if(!exists(as.character(match.call()$dat),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$dat, " not found"))
  if(!is.element(factanal, c(T, F)))  stop("Arg. must be TRUE or FALSE")
  if(method == "FA") {
    if(missing(range.Q))              stop("Arg. range.Q must be specified")
    assign("range.Q", range.Q, envir=.GlobalEnv)
  }
  if(!is.element(centering, c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(scaling, c(T, F)))   stop("Arg. must be TRUE or FALSE")
  if(!is.element(print,   c(T, F)))   stop("Arg. must be TRUE or FALSE")
  if(!is.element(profile, c(T, F)))   stop("Arg. must be TRUE or FALSE")
  assign("N", nrow(dat), envir=.GlobalEnv)
  assign("P", sum(sapply(dat, is.numeric)), envir=.GlobalEnv)
  
  # Define full conditional & Gibbs Sampler functions for desired method
  source(paste(getwd(), "/IMIFA-GIT/Gibbs_", method, ".R", sep=""))
  gibbs.arg <- list(dat, n.iters, burnin, thinning, n.store, centering, scaling, print)
  
  if(profile) {
    Rprof()
  }
    if(method == 'IFA') {
    if(missing(Q.star)) {
      Q.star       <- min(round(5 * log(P)), P)
    } else if(Q.star > P)             stop("Number of factors must be less than the number of variables")
    if(!is.element(adapt,   c(T, F))) stop("Arg. must be TRUE or FALSE")
    gibbs.arg      <- append(gibbs.arg, list(adapt))
    imifa          <- vector("list", length(Q.star))
    start.time     <- proc.time()
    imifa[[1]]     <- do.call(paste0("gibbs.", method),                          
                              args=append(list(Q.star), gibbs.arg))
  } else if(method == 'FA') {
    if((length(range.Q)  == 1 && range.Q >= P) || 
       (length(range.Q)   > 1 && any(range.Q) >= P))  
                                      stop ("Number of factors must be less than the number of variables")
    imifa          <- vector("list", length(range.Q))
    if(length(range.Q)   == 1) {
      start.time   <- proc.time()
      imifa[[1]]   <- do.call(paste0("gibbs.", method), 
                              args=append(list(range.Q), gibbs.arg))
    } else {
      for(q in range.Q) { 
        Q.ind      <- q - min(range.Q) + 1
        start.time <- proc.time()
        imifa[[Q.ind]]   <- do.call(paste0("gibbs.", method),
                                   args=append(list(q), gibbs.arg))
        cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete"))
      }
    }
  }
  tot.time  <- proc.time() - start.time
  avg.time  <- tot.time/ifelse(method == "FA", length(range.Q), length(Q.star))
  if(profile) {
    Rprof(NULL)
    print(summaryRprof())
    invisible(file.remove("Rprof.out"))
  }
    
  dat.name  <- as.character(match.call()$dat)
  attr(imifa, "Date")    <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors") <- if(method == 'FA') range.Q else Q.star
  attr(imifa, "Method")  <- method
  attr(imifa, "Name")    <- paste0(toupper(substr(dat.name, 1, 1)),
                                   substr(dat.name, 2, nchar(dat.name)))
  attr(imifa, "Time")    <- list(Total = tot.time, Average = avg.time) 
  if(print)    print(attr(imifa, "Time"))  
  attrs     <- attributes(imifa)
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac)) factanal  <- T
  if(factanal) {
    if(missing(Q.fac)) assign("Q.fac", round(sqrt(P)))
    fac     <- factanal(dat[,sapply(dat, is.numeric)], 
                         factors=Q.fac, control=list(nstart=50))
    imifa   <- append(imifa, list(fac=fac))
    attributes(imifa)    <- attrs
    names(imifa)         <- c(paste0("IMIFA", 1:(length(imifa) - 1)), "fac")
    return(imifa)
  } else {
    attributes(imifa)    <- attrs
    return(imifa)
  }
}

source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))