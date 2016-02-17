######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
 #set.seed(1)
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
  # Meat
    load(file=paste0(getwd(), "/Data/Meat.Rdata", sep=""), envir=.GlobalEnv)
    spectra  <- t(spectra); rm(last.warning)
    matplot(t(spectra), type="l", col=1:nlevels(as.factor(type)), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
  # Subjects 
    subjects <- read.csv(paste0(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  # Cereal 
    cereal   <- read.csv(paste0(getwd(), "/Data/", "Cereal.csv", sep=""))
  # Simulated data
   #SimData  <- sim.imifa(N=20, P=100, Q=15)
   #save(SimData, file=paste0(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste0(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Run the Gibbs Sampler
  sim  <- imifa(wine, method="IFA")
 #sim  <- imifa(wine, method="classify", Label=Lab)

# Save / Load results
  save(sim, file=paste0(getwd(), "/Simulations/", attr(sim, "Name"), 
                        "__Simulations_", attr(sim, "Method"), 
                        ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Simulations_", "IFA", 
                   ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  res <- tune.sims(sim)
 #res <- tune.sims(sim, Q=2)
  
# Posterior Summaries & Plots, etc.  
  plot(res, "v")
  
# Means
  plot(res, "t", "m")
  plot(res, "t", "m", mat=F)
  plot(res, "d", "m")
  plot(res, "d", "m", mat=F)
  plot(res, "p", "m")
  plot(res, "c", "m")
  
# Scores
  plot(res, "t", "s")
  plot(res, "t", "s", mat=F)
  plot(res, "d", "s")
  plot(res, "d", "s", mat=F)
  plot(res, "p", "s", Lab)
  plot(res, "c", "s")
      
# Loadings
  plot(res, "t", "l")
  plot(res, "t", "l", mat=F)
  plot(res, "d", "l")
  plot(res, "d", "l", mat=F)
  res$post.load
  plot(res, "p", "l", heat=F)
  plot(res, "p", "l")
  plot(res, "c", "l")
  
# Uniquenesses
  plot(res, "t", "u")
  plot(res, "t", "u", mat=F)
  plot(res, "d", "u")
  plot(res, "d", "u", mat=F)
  summary(res$post.psi)
  plot(res, "p", "u")
  plot(res, "c", "u")
####