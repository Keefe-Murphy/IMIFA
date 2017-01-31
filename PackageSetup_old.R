#################################################
### Set-up for Keefe Murphy's IMIFA R Package ###
#################################################

packages  <- c("abind", "corpcor", "dichromat", "e1071", "gclus", "graphics", "grDevices", "matrixcalc",
               "matrixStats", "mclust", "MCMCpack", "mvnfast", "plotrix", "Rfast", "slam", "stats", "utils")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  suppressMessages(install.packages(setdiff(packages, rownames(installed.packages()))))
}
if(length(setdiff(packages, (.packages()))) > 0) {
  suppressMessages(lapply(setdiff(packages, (.packages())), library, ch=TRUE))
}
if(!("IMIFA.env"  %in% search())) packageStartupMessage("   ________  __________________\n  /_  __/  |/   /_  __/ ___/ _ \\  \n   / / / /|_// / / / / /__/ /_\\ \\ \n _/ /_/ /   / /_/ /_/ ___/ /___\\ \\ \n/____/_/   /_/_____/_/  /_/     \\_\\    version 1.1.0")
while("IMIFA.env" %in% search())  detach("IMIFA.env")
x="..."
IMIFA.env <- new.env()
source(paste0(getwd(),   "/R/Diagnostics.R"),       local=IMIFA.env)
source(paste0(getwd(),   "/R/FullConditionals.R"),  local=IMIFA.env)
source(paste0(getwd(),   "/R/MainFunction.R"),      local=IMIFA.env)
source(paste0(getwd(),   "/R/PlottingFunctions.R"), local=IMIFA.env)
source(paste0(getwd(),   "/R/SimulateData.R"),      local=IMIFA.env)
for(meth in c("FA",   "IFA",   "MFA",  "MIFA",
              "OMFA", "OMIFA", "IMFA", "IMIFA")) {
  source(paste0(getwd(), "/R/Gibbs_", meth, ".R"),  local=IMIFA.env)
}
attach(IMIFA.env); rm(IMIFA.env, meth, packages)
