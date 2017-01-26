######################################################
### INFINITE MIXTURES OF INFINITE FACTOR ANALYZERS ###
######################################################

# Preamble
 #set.seed(1)
 #rm(list=ls(all=TRUE))
 #while("IMIFA.env" %in% search()) detach("IMIFA.env")
  if(!is.element(getwd(), c("/home/kmurphy", "/home/kmurphy/IMIFA"))) {
    wd       <- try(setwd("C:/Users/Windows/Documents/IMIFA/"), silent=TRUE)
    if(inherits(wd, "try-error")) {
      setwd("D:/Documents/IMIFA"); rm(wd)
      datdir <- "D:/Dropbox/UCD/IMIFA"
    } else {
      datdir <- "C:/Users/Windows/Dropbox/UCD/IMIFA"
    }
  } else   {
    datdir   <- "/home/kmurphy"
    setwd("/home/kmurphy/IMIFA-GIT")
  }
  source(paste0(getwd(), "/PackageSetup_old.R"))

# Read in the data
  # Olive       (# area, region)
    load(file=paste0(datdir, "/Data/Olive.Rdata", sep=""), envir=.GlobalEnv)
    area     <- as.factor(olive$area)
    cinzia   <- ifelse(olive$region < 5, 1, ifelse(olive$region < 7, 2, ifelse(olive$region == 9, 4, 3)))
    region   <- as.factor(olive$region)
    olive    <- olive[,-c(1, 2)]
    # Robust Olive
      load(file=paste0(datdir, "/Data/ExtraOlive.Rdata", sep=""), envir=.GlobalEnv)
  # Wine        (# lab)
    load(file=paste0(datdir, "/Data/Wine.Rdata", sep=""), envir=.GlobalEnv)
    lab      <- wine[,1]
    # 13 Variable Version [% wine13]
      load(file=paste0(datdir, "/Data/Wine13.Rdata", sep=""), envir=.GlobalEnv)
      lab    <- wine13[,1]
  # Iris        (# species)
    load(file=paste0(datdir, "/Data/Iris.Rdata", sep=""), envir=.GlobalEnv)
    species  <- iris[,5]
  # Oils        (# oliveoillabels)
    load(file=paste0(datdir, "/Data/Oliveoils.Rdata", sep=""), envir=.GlobalEnv)
    oils     <- t(oliveoils); rm(oliveoils, wavelengths, x)
  # Coffee      (# type, country) [% coffee]
    load(file=paste0(datdir, "/Data/Coffee.Rdata", sep=""), envir=.GlobalEnv)
  # Brain       (# region) {n.b. pareto scaling}
    load(file=paste0(datdir, "/Data/RatBrain.Rdata", sep=""), envir=.GlobalEnv)
   #library(MetabolAnalyze); data(BrainSpectra)
    region   <- brain$Region;      brain  <- brain[,-1]
   #brain    <- BrainSpectra[[1]]; region <- BrainSpectra[[2]]
    ppm.g    <- do.call(cbind, by(brain, region, colMeans))
   #graphics::matplot(ppm.g, type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1)
    graphics::matplot(t(brain), type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1, col=region)
    graphics::axis(1, at=seq(from=0, to=nrow(ppm.g), by=20), labels=9:1, tick=TRUE, lwd.ticks=1, xpd=TRUE)
    graphics::axis(1, at=seq_len(nrow(ppm.g)), labels=FALSE, tick=TRUE, tcl=-0.2)
    graphics::legend("topleft", cex=0.8, legend=c("Brain Stem", "Cerebellum", "Hippocampus", "Pre-frontal Cores"), bty="n", lty=1, col=seq_len(4))
  # Urine       (# grp) {n.b. pareto scaling}
    load(file=paste0(datdir, "/Data/Epi_urine_data.Rdata", sep=""), envir=.GlobalEnv)
    rat.cov  <- read.csv(file=paste0(datdir, "/Data/D10Weight.csv"))
    ppms     <- substr(colnames(x10[,4:ncol(x10)]), 2, 6); rm(x)
    urine    <- x10[,4:ncol(x10)]
    grp      <- x10[,"Group"]
    ppm.g    <- do.call(cbind, by(urine, grp, colMeans))
   #graphics::matplot(ppm.g, type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1)
    graphics::matplot(t(urine), type="l", xlab="Chemical Shift (ppm)", yaxt="n", ylab="", bty="n", xaxt="n", lwd=2, lty=1, col=grp)
    graphics::axis(1, at=seq(from=20, to=nrow(ppm.g), by=20), labels=9:1, tick=TRUE, lwd.ticks=1, xpd=TRUE)
    graphics::axis(1, at=seq_len(nrow(ppm.g)), labels=FALSE, tick=TRUE, tcl=-0.2)
    graphics::legend("topleft", legend=c("Control", "Epileptic"), bty="n", lty=1, col=c(1, 2))
    rat.cov  <- read.csv(file=paste0(datdir, "/Data/D10Weight.csv"))
    graphics::plot(rat.cov$Weight, type="h", main="", ylab="Weight (g)", xlab="Observation", xaxt="n", lty=1, lwd=2, col=replace(grp, 13, 0) + 2)
    graphics::axis(1, at=seq_len(nrow(urine)), labels=seq_len(nrow(urine)), tick=FALSE)
    graphics::legend("topright", legend=c("Control", "Epileptic", "Outlier"), bty="n", lty=1, lwd=2, col=c(3, 4, 2))
  # Meat        (# type)
    load(file=paste0(datdir, "/Data/Meat.Rdata", sep=""), envir=.GlobalEnv)
    spectra  <- t(spectra); rm(last.warning)
    meat.col <- as.numeric(as.factor(type))
    # All Meats
      graphics::matplot(t(spectra), type="l", col=meat.col, xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      graphics::legend("topleft", legend=levels(factor(type))[Rfast::sort_unique(meat.col)], bty="n", lty=1, col=sort(unique(meat.col)))
    # Red vs. White
      meat.col[meat.col == 3] <- 1; meat.col[meat.col > 1] <- 2
      graphics::matplot(t(spectra), type="l", col=meat.col, xlab="Wavelength", ylab="Spectral Reflectance", main="Meat Data")
      graphics::legend("topleft", legend=c("White Meat", "Red Meat"), bty="n", lty=1, col=c(2, 1))
  # Microarray  (# labels)
    load(file=paste0(datdir, "/Data/Microarray.Rdata", sep=""), envir=.GlobalEnv)
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
    subjects <- read.csv(paste0(datdir, "/Data/", "SubjectMarks.csv", sep=""))
  # Cereal      (# classes)
    cereal   <- read.csv(paste0(datdir, "/Data/", "Cereal.csv", sep=""))
    classes  <- cereal$Cereals
  # Simulated   (# classes)
   #SimData  <- sim_IMIFA_data(N=80, G=4, P=100, Q=c(5, 1, 4, 0), nn=c(20, 20, 20, 20))
   #save(SimData, file=paste0(datdir,"/Data/Simulated_Data.Rdata", sep=""))
    load(file=paste0(datdir, "/Data/Simulated Data/Replication 1/Simulated_DataN25P50.Rdata", sep=""), envir=.GlobalEnv)
    classes  <- attr(SimData, "Labels")

# Run the Gibbs Sampler
  sim        <- mcmc_IMIFA(olive, method="IMIFA")
  summary(sim)

# Save / Load Simulations
  save(sim, file=paste0(datdir, "/Simulations/", attr(sim, "Name"),
                        "__Simulations_", attr(sim, "Method"), ".Rdata"))
  load(file=paste0(datdir, "/Simulations/", "wine",
                   "__Simulations_", "IMIFA", ".Rdata"), envir=.GlobalEnv)

# Posterior Summaries (optional: additional 'burnin' & 'thinning', user-defined G/Q, model selection criterion)
  res        <- get_IMIFA_results(sim, zlabels=area)
  summary(res)

# Save / Load Results
  save(res, file=paste0(datdir, "/Simulations/", attr(res, "Name"),
                        "__Results_", attr(res, "Method"), ".Rdata"))
  load(file=paste0(datdir, "/Simulations/", "wine",
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
