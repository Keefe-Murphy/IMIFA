######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
 #set.seed(1)
 #rm(list=ls(all=TRUE))
 #while("IMIFA.env" %in% search()) detach("IMIFA.env")
  if(getwd() != "/home/kmurphy")  {
    wd       <- try(setwd("C:/Users/Windows/Dropbox/UCD/Claire IMIFA"), silent=TRUE)
    if(inherits(wd, "try-error")) setwd("D:/Dropbox/UCD/Claire IMIFA"); rm(wd)
  }
  source(paste0(getwd(), "/IMIFA-GIT/PackageSetup.R"))
    
# Read in the data
  # Wine        (# lab)
    load(file=paste0(getwd(), "/Data/Wine.Rdata", sep=""), envir=.GlobalEnv)
    lab      <- wine[,1]
    # 13 Variable Version [% wine13]
      load(file=paste0(getwd(), "/Data/Wine13.Rdata", sep=""), envir=.GlobalEnv)
      lab    <- wine13[,1]
  # Iris        (# species)
    load(file=paste0(getwd(), "/Data/Iris.Rdata", sep=""), envir=.GlobalEnv)
    species  <- iris[,5]
  # Olive       (# area, region)
    load(file=paste0(getwd(), "/Data/Olive.Rdata", sep=""), envir=.GlobalEnv)
    area     <- as.factor(olive$area)
    region   <- as.factor(olive$region)
    olive    <- olive[,-c(1, 2)]
    # Oils      (# oliveoillabels)
      load(file=paste0(getwd(), "/Data/Oliveoils.Rdata", sep=""), envir=.GlobalEnv)
      oils   <- t(oliveoils); rm(oliveoils, wavelengths, x)
  # Coffee      (# type, country) [% coffee]
    load(file=paste0(getwd(), "/Data/Coffee.Rdata", sep=""), envir=.GlobalEnv) 
  # Urine       (# grp) {n.b. pareto scaling}
    load(file=paste0(getwd(), "/Data/Epi_urine_data.Rdata", sep=""), envir=.GlobalEnv)
    ppms     <- substr(colnames(x10[,4:ncol(x10)]), 2, 6); rm(x)
    urine    <- x10[,4:ncol(x10)]
    grp      <- x10[,"Group"]
    ppm.g    <- vapply(seq_len(max(grp)), function(g) colMeans(urine[grp == g,]), numeric(ncol(urine)))
    matplot(ppm.g, type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1)
    axis(1, at=seq(from=20, to=nrow(ppm.g), by=20), labels=9:1, tick=TRUE, lwd.ticks=1, xpd=TRUE)
    axis(1, at=seq_len(nrow(ppm.g)), labels=FALSE, tick=TRUE, tcl=-0.2)
    legend("topleft", legend=c("Control", "Epileptic"), bty="n", lty=1, col=c(1, 2))  
  # Meat        (# type)
    load(file=paste0(getwd(), "/Data/Meat.Rdata", sep=""), envir=.GlobalEnv)
    spectra  <- t(spectra); rm(last.warning)
    # All Meats
      matplot(t(spectra), type="l", col=c(2, 1, 3, 4, 5), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      legend("topleft", legend=levels(factor(type)), bty="n", lty=1, col=c(2, 1, 3, 4, 5))
    # Red vs. White
      matplot(t(spectra), type="l", col=c(2, 1, 2, 1, 1), xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      legend("topleft", legend=c("Red Meat", "White Meat"), bty="n", lty=1, col=c(2, 1))
  # Microarray  (# labels)
    load(file=paste0(getwd(), "/Data/Microarray.Rdata", sep=""), envir=.GlobalEnv)
    aliza    <- data.frame(labels = as.factor(mydata.alizadeth$y), t(mydata.alizadeth$x))
    alon     <- data.frame(labels = as.factor(mydata.alon$y),        mydata.alon$x)
    golub    <- data.frame(labels = as.factor(mydata.golub$y),     t(mydata.golub$x))
    khan     <- data.frame(labels = as.factor(mydata.khan$y),      t(mydata.khan$x))
    mammary  <- data.frame(labels = as.factor(mydata.mammary$y),   t(mydata.mammary$x))
    # Remove 'mydata.' objects
      remove <- ls()
      remove <- c(remove[grepl("^mydata.", remove)], "remove")
     #remove <- c(remove, setdiff(c("aliza", "alon", "golub", "khan", "mammary"), "golub"))
      rm(list = remove)
      labels <- as.factor(as.matrix(golub[1]))
     #G2levs <- ifelse(labels == "AML", "AML", "ALL")
  # Subjects 
    subjects <- read.csv(paste0(getwd(), "/Data/", "SubjectMarks.csv", sep=""))
  # Cereal      (# classes)
    cereal   <- read.csv(paste0(getwd(), "/Data/", "Cereal.csv", sep=""))
    classes  <- cereal$Cereals
  # Simulated   (# classes)
   #SimData  <- sim.IMIFA(N=80, G=4, P=100, Q=c(5, 1, 4, 0), nn=c(20, 20, 20, 20))
   #save(SimData, file=paste0(getwd(),"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste0(getwd(), "/Data/Simulated_Data.Rdata", sep=""), envir=.GlobalEnv)
    classes  <- attr(SimData, "Labels")

# Run the Gibbs Sampler
  sim        <- mcmc.IMIFA(wine, method="IMIFA")
  summary(sim)

# Save / Load Simulations
  save(sim, file=paste0(getwd(), "/Simulations/", attr(sim, "Name"), 
                        "__Simulations_", attr(sim, "Method"), ".Rdata"))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Simulations_", "IMIFA", ".Rdata"), envir=.GlobalEnv)
                        
# Posterior Summaries (optional: additional 'burnin' & 'thinning', user-defined G/Q, model selection criterion)
  res        <- tune.IMIFA(sim, zlabels=lab)
  summary(res)

# Save / Load Results
  save(res, file=paste0(getwd(), "/Simulations/", attr(res, "Name"), 
                        "__Results_", attr(res, "Method"), ".Rdata"))
  load(file=paste0(getwd(), "/Simulations/", "wine", 
                   "__Results_", "IMIFA", ".Rdata"), envir=.GlobalEnv)

# Model Selection Parameters
  plot(res, "GQ") 

# Cluster Labels
  plot(res, "z")
  res$Clust$map
  res$Clust$lab.rate

# Means
  plot(res, "p", "m")
  plot(res, "a", "m")
  plot(res, "a", "m", mat=FALSE)
  plot(res, "t", "m")
  plot(res, "t", "m", mat=FALSE)
  plot(res, "d", "m")
  plot(res, "d", "m", mat=FALSE)
  plot(res, "m", "m")
  plot(res, "c", "m")
  res$Means$post.mu
  res$Means$var.mu
  res$Means$ci.mu
  
# Loadings
  plot(res, "p", "l")
  plot(res, "a", "l")
  plot(res, "a", "l", mat=FALSE)
  plot(res, "a", "l", load.meth="raw")
  plot(res, "t", "l")
  plot(res, "t", "l", mat=FALSE)
  plot(res, "d", "l")
  plot(res, "d", "l", mat=FALSE)
  plot(res, "m", "l")
  plot(res, "m", "l", load.meth="raw")
  plot(res, "c", "l")
  res$Loadings$post.load
  res$Loadings$var.load
  res$Loadings$ci.load

# Scores
  plot(res, "a", "s")
  plot(res, "a", "s", lab)
  plot(res, "a", "s", mat=FALSE)
  plot(res, "t", "s")
  plot(res, "t", "s", mat=FALSE)
  plot(res, "d", "s")
  plot(res, "d", "s", mat=FALSE)
  plot(res, "m", "s")
  plot(res, "c", "s")
  res$Scores$post.eta
  res$Scores$var.eta
  res$Scores$ci.eta
  
# Uniquenesses
  plot(res, "p", "u")
  plot(res, "a", "u")
  plot(res, "a", "u", mat=FALSE)
  plot(res, "t", "u")
  plot(res, "t", "u", mat=FALSE)
  plot(res, "d", "u")
  plot(res, "d", "u", mat=FALSE)
  plot(res, "m", "u")
  plot(res, "c", "u")
  res$Uniquenesses$post.psi
  res$Uniquenesses$var.psi
  res$Uniquenesses$ci.psi

# Mixing Proportions
  plot(res, "a", "p")
  plot(res, "a", "p", mat=FALSE)
  plot(res, "t", "p")
  plot(res, "t", "p", mat=FALSE)
  plot(res, "d", "p")
  plot(res, "d", "p", mat=FALSE)
  plot(res, "m", "p")
  plot(res, "m", "p", mat=FALSE)
  plot(res, "c", "p")
  res$Clust$post.pi
  res$Clust$var.pi
  res$Clust$ci.pi
  
# Covariance Matrices & Error Metrics
  plot(res, "e")
  
# Learning 'alpha' for Dirichlet Process methods
  plot(res, "a", "a")
  plot(res, "t", "a")
  plot(res, "d", "a")
  plot(res, "m", "a")
  plot(res, "c", "a")
#####