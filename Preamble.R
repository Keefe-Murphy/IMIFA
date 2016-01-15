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

initialise <- function(data, factanal=T, Q=NULL, range.Q=NULL, Q.fac=NULL,
                       method=c("IMIFA", "MIFA", "IFA", "FA"), ...) {
  
  method   <- match.arg(method)
  assign("method", method, envir=.GlobalEnv)
  if(!is.element(factanal, c(T, F)))  stop("Arg. must be TRUE or FALSE")
  if(method == "FA" && missing(range.Q)) stop("Arg. range.Q must be specified")
  assign("range.Q", range.Q, envir=.GlobalEnv)
  
# Define full conditional & Gibbs Sampler functions for desired method
  source(paste(dataDirectory, "/IMIFA-GIT/Gibbs_", method, ".R", sep=""))

# Vanilla 'factanal' for comparison purposes
  if(factanal) {
    if(missing(Q.fac)) assign("Q.fac", round(sqrt(sum(sapply(data, is.numeric)))))
    if(!exists("Label")) stop("Should the data be labelled?")
    fac   <- factanal(data[,sapply(data, is.numeric)], 
                         factors=Q.fac, control=list(nstart=50))
  }
  return(list(fac=fac))
}