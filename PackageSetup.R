#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages  <- c("MCMCpack", "slam")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  invisible(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  invisible(lapply(setdiff(packages, (.packages())), library, ch=T))
}
rm(packages)

imifa.gibbs <- function(dat=NULL, method=c("IMIFA", "MIFA", "MFA", "IFA", "FA", "classify"), n.iters=50000,
                        Label=NULL, factanal=F, Q.star=NULL, range.Q=NULL, Q.fac=NULL, thinning=2,
                        burnin=n.iters/5 - 1, n.store=ceiling((n.iters - burnin)/thinning),
                        centering=T, scaling=T, print=T, adapt=T,
                        sigma.mu=NULL, sigma.l=NULL, psi.alpha=NULL, psi.beta=NULL,
                        phi.nu=NULL, delta.a1=NULL, delta.a2=NULL, profile=F, ...) {
  
  method    <- match.arg(method)
  if(missing(dat))                    stop("Dataset must be supplied")
  if(!exists(as.character(match.call()$dat),
             envir=.GlobalEnv))       stop(paste0("Object ", match.call()$dat, " not found"))
  if(!is.element(factanal,  c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(centering, c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(scaling,   c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(print,     c(T, F))) stop("Arg. must be TRUE or FALSE")
  if(!is.element(profile,   c(T, F))) stop("Arg. must be TRUE or FALSE")
  
  # Remove non-numeric columns & apply centering & scaling if necessary 
  dat       <- as.data.frame(dat)
  dat       <- dat[sapply(dat, is.numeric)]
  dat       <- scale(dat, center=centering, scale=scaling)
  
  # Define full conditionals, hyperparamters & Gibbs Sampler function for desired method
  N         <- nrow(dat)
  P         <- ncol(dat)
  if(is.null(rownames(dat))) rownames(dat) <- c(1:N)
  if(missing("sigma.mu"))    sigma.mu      <- 0.5
  if(missing("psi.alpha"))   psi.alpha     <- 2
  if(missing("psi.beta"))    psi.beta      <- 0.6
  if(method == "FA") {
    if(missing("sigma.l"))   sigma.l       <- 0.5
  } else if(method == "IFA" ||
            method == "classify") {
    if(missing("phi.nu"))    phi.nu        <- 3
    if(missing("delta.a1"))  delta.a1      <- 2.1
    if(missing("delta.a2"))  delta.a2      <- 12.1
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
                    thinning=thinning, n.store=n.store, print=print)
  if(profile)  Rprof()
  if(method == "IFA" ||
     method == "classify") {
     gibbs.arg     <- append(gibbs.arg, list(phi.nu=phi.nu, delta.a1=delta.a1, delta.a2=delta.a2,
                                             adapt=adapt, b0=0.1, b1=0.00005, prop=3/4, 
                                             epsilon=ifelse(centering, 0.1, 0.01)))
  if(missing(Q.star)) {
     Q.star        <- min(floor(3 * log(P)), P)
  } else if(Q.star  > P)              stop("Number of factors must be less than the number of variables")
    if(!is.element(adapt,   c(T, F))) stop("Arg. must be TRUE or FALSE")
    if(method == "IFA") {
      imifa        <- vector("list", length(Q.star))
      start.time   <- proc.time()
      imifa[[1]]   <- do.call(paste0("gibbs.", method),                          
                              args=append(list(data=dat, N=N, Q=Q.star), gibbs.arg))
    } else {
      if(missing(Label))              stop("Data must be labelled for classification")
      if(!exists(as.character(match.call()$Label),
                 envir=.GlobalEnv))   stop(paste0("Object ", match.call()$Label, " not found"))
      Label   <- as.factor(Label)
      if(length(Label) != N)          stop(paste0("Labels must be a factor of length N=",  n.obs))
      imifa   <- vector("list", nlevels(Label))
      start.time   <- proc.time()
      for(i in 1:nlevels(Label)) {
        temp.dat   <- dat[Label==levels(Label)[i],]
        imifa[[i]] <- do.call(paste0("gibbs.", "IFA"),
                              args=append(list(data=temp.dat, N=nrow(temp.dat), Q=Q.star), gibbs.arg))
        cat(paste0(round(i/nlevels(Label) * 100, 2), "% Complete\n"))
      }
    }
  } else if(method == "FA") {
    gibbs.arg      <- append(gibbs.arg, list(sigma.l))
    if(missing(range.Q))              stop("Arg. range.Q must be specified")
    if((length(range.Q)  == 1 && range.Q >= P) || 
       (length(range.Q)   > 1 && any(range.Q) >= P))  
                                      stop ("Number of factors must be less than the number of variables")
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
        cat(paste0(round(Q.ind/length(range.Q) * 100, 2), "% Complete\n"))
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
  attr(imifa, "Date")    <- format(Sys.Date(), "%d-%b-%Y")
  attr(imifa, "Factors") <- if(method == "FA") range.Q else Q.star
  attr(imifa, "Method")  <- paste0(toupper(substr(method, 1, 1)),
                                   substr(method, 2, nchar(method)))
  attr(imifa, "Name")    <- paste0(toupper(substr(dat.name, 1, 1)),
                                   substr(dat.name, 2, nchar(dat.name)))
  attr(imifa, "Store")   <- n.store
  attr(imifa, "Time")    <- list(Total = tot.time, Average = avg.time) 
  if(print)    print(attr(imifa, "Time"))  
      
# Vanilla 'factanal' for comparison purposes
  if(!missing(Q.fac)) factanal  <- T
  if(factanal) {
    if(missing(Q.fac)) Q.fac    <- round(sqrt(P))
    fac     <- factanal(dat, factors=Q.fac, control=list(nstart=50))
    imifa   <- append(imifa, list(fac=fac))
    names(imifa)         <- c(paste0("IMIFA", 1:(length(imifa) - 1)), "fac")
  } 
  class(imifa)           <- "IMIFA"
  return(imifa)
}

source(paste(getwd(), "/IMIFA-GIT/Diagnostics.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/PlottingFunctions.R", sep=""))
source(paste(getwd(), "/IMIFA-GIT/Simulate_Data.R", sep=""))