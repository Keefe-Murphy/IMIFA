######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
 #set.seed(1)
 #rm(list=ls(all=T))
  if(getwd() != "/home/kmurphy")  {
    wd       <- try(setwd("C:/Users/Windows/Dropbox/UCD/Claire IMIFA"), silent=T)
    if(inherits(wd, "try-error")) {
      setwd("D:/Dropbox/UCD/Claire IMIFA")
    }
    rm(wd)
  }
  source(paste0(getwd(), "/IMIFA-GIT/PackageSetup.R", sep=""))
    
# Read in the data
  # Wine
    load(file=paste0(getwd(), "/Data/Wine.Rdata", sep=""), envir=.GlobalEnv)
    Lab      <- wine[,1]
  # Iris
    load(file=paste0(getwd(), "/Data/Iris.Rdata", sep=""), envir=.GlobalEnv)
    Species  <- iris[,5]
  # Urine
    load(file=paste0(getwd(), "/Data/Epi_urine_data.Rdata", sep=""), envir=.GlobalEnv)
    ppms     <- substr(colnames(x10[,4:ncol(x10)]), 2,6); rm(x)
    urine    <- x10[,4:ncol(x10)]
    Grp      <- x10[,"Group"]
    ppm.g    <- do.call(cbind, lapply(seq_len(max(Grp)), function(g) colMeans(urine[Grp == g,])))
    matplot(ppm.g, type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1)
    axis(1, at=seq(from=20, to=nrow(ppm.g), by=20), labels=9:1, tick=T, lwd.ticks=1, xpd=T)
    axis(1, at=seq_len(nrow(ppm.g)), labels=FALSE, tick=T, tcl=-0.2)
    legend("topleft", legend=c("Control", "Epileptic"), bty="n", lty=1, col=c(1,2))  
  # Meat
    load(file=paste0(getwd(), "/Data/Meat.Rdata", sep=""), envir=.GlobalEnv)
    spectra  <- t(spectra); rm(last.warning)
    # All Meats
      matplot(t(spectra), type="l", col=c(2,1,3,4,5), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      legend("topleft", legend=levels(factor(type)), bty="n", lty=1, col=c(2,1,3,4,5))
    # Red vs. White
      matplot(t(spectra), type="l", col=c(2,1,2,1,1), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      legend("topleft", legend=c("Red Meat", "White Meat"), bty="n", lty=1, col=c(2,1))
  # Microarray
    load(file=paste0(getwd(), "/Data/Microarray.Rdata", sep=""), envir=.GlobalEnv)
    aliza    <- data.frame(as.factor(mydata.alizadeth$y), t(mydata.alizadeth$x))
    golub    <- data.frame(as.factor(mydata.golub$y), t(mydata.golub$x))
    khan     <- data.frame(as.factor(mydata.khan$y), t(mydata.khan$x))
    mammary  <- data.frame(as.factor(mydata.mammary$y), t(mydata.mammary$x))
    # Remove 'mydata.' objects
      remove <- ls()
      remove <- c(remove[grepl("^mydata.", remove)], "remove")
     #remove <- c(remove, setdiff(c("aliza", "golub", "khan", "mammary"), "aliza"))
      rm(list=remove)
  # Subjects 
    subjects <- read.csv(paste0(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  # Cereal 
    cereal   <- read.csv(paste0(getwd(), "/Data/", "Cereal.csv", sep=""))
  # Simulated data
   #SimData  <- sim.IMIFA(N=80, G=4, P=100, Q=c(5, 1, 4, 0), nn=c(20, 20, 20, 20))
   #save(SimData, file=paste0(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste0(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)
    cl <- attr(SimData, "Labels")

# Run the Gibbs Sampler
  sim  <- IMIFA.mcmc(wine, method="IMFA", range.Q=3)
 #sim  <- IMIFA.mcmc(wine, method="classify", Label=Lab)

# Save / Load Simulations
  save(sim, file=paste0(getwd(), "/Simulations/", attr(sim, "Name"), 
                        "__Simulations_", attr(sim, "Method"), 
                        ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Simulations_", "OMIFA", 
                   ".Rdata", sep=""), envir=.GlobalEnv)

# Posterior Summaries (optional: additional 'burnin' & 'thinning', user-defined G/Q, model selection criterion)
  res <- tune.IMIFA(sim, Labels=Lab)
 #res <- tune.IMIFA(sim, G=3, Q=3, criterion="bicm", Labels=Lab)
  res$Error

# Save / Load Results
  save(res, file=paste0(getwd(), "/Simulations/", attr(res, "Name"), 
                        "__Results_", attr(res, "Method"), 
                        ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Results_", "OMIFA", 
                   ".Rdata", sep=""), envir=.GlobalEnv)

# Model Selection Parameters
  plot(res, "GQ") 

# Cluster Labels
  plot(res, "Z")

# Means
  plot(res, "a", "m")
  plot(res, "a", "m", mat=F)
  plot(res, "t", "m")
  plot(res, "t", "m", mat=F)
  plot(res, "d", "m")
  plot(res, "d", "m", mat=F)
  do.call(rbind, lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$post.mu))
  do.call(rbind, lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$var.mu))
  plot(res, "p", "m")
  plot(res, "c", "m")
  
# Scores
  plot(res, "a", "s")
  plot(res, "a", "s", Lab)
  plot(res, "a", "s", mat=F)
  plot(res, "t", "s")
  plot(res, "t", "s", mat=F)
  plot(res, "d", "s")
  plot(res, "d", "s", mat=F)
  plot(res, "p", "s")
  plot(res, "c", "s")
      
# Loadings
  plot(res, "a", "l")
  plot(res, "a", "l", mat=F)
  plot(res, "t", "l")
  plot(res, "t", "l", mat=F)
  plot(res, "d", "l")
  plot(res, "d", "l", mat=F)
  lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$post.load)
  lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$var.load)
  plot(res, "p", "l")
  plot(res, "c", "l")
  
# Uniquenesses
  plot(res, "a", "u")
  plot(res, "a", "u", mat=F)
  plot(res, "t", "u")
  plot(res, "t", "u", mat=F)
  plot(res, "d", "u")
  plot(res, "d", "u", mat=F)
  do.call(rbind, lapply(seq_len(res$GQ.results$G), function(g) summary(res[[g]]$post.psi)))
  do.call(rbind, lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$var.psi))
  plot(res, "p", "u")
  plot(res, "c", "u")

# Mixing Proportions
  
# Covariance Matrices
  lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$cov.est)

# Metropolis Hastings 'alpha' for Dirichlet Process methods
  plot(res, "a", "a")
  plot(res, "t", "a")
  plot(res, "d", "a")
  plot(res, "p", "a")
  plot(res, "c", "a")
####