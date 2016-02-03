######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
  if(getwd() != "/home/kmurphy") {
    setwd("C:/Users/Windows/Documents/Claire IMIFA")
  }
  source(paste0(getwd(), "/IMIFA-GIT/PackageSetup.R", sep=""))
    
# Read in the data
  load(file=paste0(getwd(), "/Data/Wine.Rdata", sep=""), envir=.GlobalEnv)
  Lab        <- wine[,1]
  #load(file=paste0(getwd(), "/Data/Meat.Rdata", sep=""), envir=.GlobalEnv)
  #spectra   <- t(spectra); rm(last.warning)
  #matplot(t(spectra), type="l", col=1:nlevels(as.factor(type)), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
  #subjects  <- read.csv(paste0(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  #cereal    <- read.csv(paste0(getwd(), "/Data/", "Cereal.csv", sep=""))

# Simulate data
  SimData    <- sim.imifa()
  #save(SimData, file=paste0(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
  load(file=paste0(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)

# Run the Gibbs Sampler
  sim  <- imifa.gibbs(wine, method="IFA")
  #sim <- imifa.gibbs(wine, method="classify", Label=Lab)

# Save / Load results
  save(sim, file=paste0(getwd(), "/Simulations/", 
                        attr(sim, "Name"), "__Simulations_", attr(sim, "Method"), "_", 
                        attr(sim, "Date"), ".Rdata", sep=""))
  load(file=paste0(getwd(), "/Simulations/", 
                   attr(sim, "Name"), "__Simulations_", attr(sim, "Method"), "_", 
                   attr(sim, "Date"), ".Rdata", sep=""), envir=.GlobalEnv)

# Convergence diagnostics (optional: additional 'burnin' & 'thinning' & user-defined Q)
  res <- tune.sims(sim)
 #res <- tune.sims(sim, Q=2)
  
# Posterior Summaries & Plots, etc.  
  plot(res, "c")
  
# Means
  plot(res, "t", "m")
  plot(res, "t", "m", mat=F)
  plot(res, "p", "m")
  plot(res, "a", "m")
  
# Scores
  plot(res, "t", "s")
  plot(res, "t", "s", mat=F)
  plot(res, "p", "s", Lab)
  plot(res, "a", "s")
      
# Loadings
  plot(res, "t", "l")
  plot(res, "t", "l", mat=F)
  plot(res, "p", "l", heat=F)
  plot(res, "p", "l")
  plot(res, "a", "l")
  
# Uniquenesses
  plot(res, "t", "u")
  plot(res, "t", "u", mat=F)
  plot(res, "p", "u")
  plot(res, "a", "u")
####