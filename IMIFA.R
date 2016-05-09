######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
 #set.seed(1)
 #rm(list=ls(all=T))
  if(getwd() != "/home/kmurphy") {
    setwd("C:/Users/Windows/Documents/Claire IMIFA")
  }
  source(paste0(getwd(), "/IMIFA-GIT/PackageSetup.R", sep=""))
    
# Read in the data
  # Wine
    load(file=paste0(getwd(), "/Data/Wine.Rdata", sep=""), envir=.GlobalEnv)
    Lab      <- wine[,1]
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
  # Subjects 
    subjects <- read.csv(paste0(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  # Cereal 
    cereal   <- read.csv(paste0(getwd(), "/Data/", "Cereal.csv", sep=""))
  # Simulated data
   #SimData  <- sim.imifa(N=80, G=4, P=100, Q=c(5, 1, 4, 0), nn=c(20, 20, 20, 20))
   #save(SimData, file=paste0(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste0(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Run the Gibbs Sampler
  sim  <- imifa.mcmc(wine, method="MIFA", range.G=3, z.list=Lab)
 #sim  <- imifa.mcmc(wine, method="classify", Label=Lab)

# Save / Load Simulations
  save(sim, file=paste0(getwd(), "/Simulations/", attr(sim, "Name"), 
                        "__Simulations_", attr(sim, "Method"), 
                        ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Simulations_", "MIFA", 
                   ".Rdata", sep=""), envir=.GlobalEnv)

# Posterior Summaries (optional: additional 'burnin' & 'thinning', user-defined G/Q, model selection criterion)
  res <- tune.imifa(sim, Labels=Lab)
  res$Clust$conf.mat
  res$Error
  plot(res, "GQ")
 #res <- tune.imifa(sim, G=3, Q=3, criterion="bicm", Labels=Lab)
  
# Save / Load Results
  save(res, file=paste0(getwd(), "/Simulations/", attr(res, "Name"), 
                        "__Results_", attr(res, "Method"), 
                        ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Results_", "MIFA", 
                   ".Rdata", sep=""), envir=.GlobalEnv)

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

# Covariance Matrices
  lapply(seq_len(res$GQ.results$G), function(g) res[[g]]$cov.est)

# Mixing Proportions

# Cluster Labels
####