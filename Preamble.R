preamble   <- function(seed=21092015, rem.lib=F, rem.all=F, ...) {
  
  if(!is.element(rem.lib, c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(rem.all, c(T, F))) stop("Arg. must be TRUE or FALSE")
  set.seed(seed)
  assign("def.par", par(), envir=.GlobalEnv)
  packages <- c("pgmm", "car", "MCMCpack", "compiler")
  if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
    invisible(install.packages(setdiff(packages, rownames(installed.packages()))))
  }
  if(length(setdiff(packages, .packages(all.available=T))) > 0) {
    invisible(lapply(setdiff(packages, .packages(all.available=T)), library, ch=T))
  }
  
# WARNING: Remove loaded libraries
  if(rem.lib) {
    pkgs   <- names(sessionInfo()$otherPkgs)
    pkgs   <- paste('package:', pkgs, sep = "")
    invisible(lapply(pkgs, detach, ch = T, unload = T, force= T))
  }

# WARNING: Remove everything
  if(rem.all) rm(list = ls(all = TRUE))
}

preamble()

imifa.gibbs <- function(dat=NULL, method=c("IMIFA", "MIFA", "IFA", "FA"), 
                       factanal=F, Q.star=NULL, range.Q=NULL, Q.fac=NULL,
                       sigma.mu=NULL, sigma.l=NULL, psi.alpha=NULL, psi.beta=NULL,
                       phi.nu=NULL, delta.a1=NULL, delta.a2=NULL, profile=F, ...) {
  
  method   <- match.arg(method)
  assign("method", method, envir=.GlobalEnv)
  if(missing(dat))                    stop("Dataset must be supplied")
  if(!is.element(factanal, c(T, F)))  stop("Arg. must be TRUE or FALSE")
  if(method == "FA") {
    if(missing(range.Q))              stop("Arg. range.Q must be specified")
    assign("range.Q", range.Q, envir=.GlobalEnv)
  }
  assign("N", nrow(dat), envir=.GlobalEnv)
  assign("P", sum(sapply(dat, is.numeric)), envir=.GlobalEnv)
  
# Define full conditional & Gibbs Sampler functions for desired method
  source(paste(dataDirectory, "/IMIFA-GIT/Gibbs_", method, ".R", sep=""))

# Vanilla 'factanal' for comparison purposes
  if(factanal) {
    if(missing(Q.fac)) assign("Q.fac", round(sqrt(P)))
    if(!exists("Label")) stop("Should the data be labelled?")
    fac   <- factanal(dat[,sapply(dat, is.numeric)], 
                         factors=Q.fac, control=list(nstart=50))
    return(list(fac=fac))
  }
}

lappend    <- function(...) {
  lists    <- list(...)
  if(length(as.list(match.call())) < 1)                  stop("lappend needs at least one list as inputs")
  if(all(lapply(lists, function(x) class(x)) != "list")) stop("At least one argument must be a list") 
  n        <- unique(unlist(lapply(lists, names)))
  names(n) <- n
  lapply(n, function(ni) unlist(lapply(lists, `[[`, ni)))
}