# IMIFA v1.1.1 - (__2nd release: 2017-02-15__)
============================================

## New Features
* Sped up simulating cluster labels from unnormalised log probabilities using the Gumbel-Max trick (Yellott, 1977).  
  \cr As a result `gumbel_max` replaces earlier function to sample cluster labels & is now unhidden/exported w/ accompanying documentation.
* Simulated finite/overfitted mixing proportions using `rgamma(G, alpha, 1)` trick (Devroye 1986, p.594) instead of `MCMCpack:rdirichlet`.  
  \cr As a result `rDirichlet` replaces earlier function to sample mixing proportions & is now unhidden/exported w/ accompanying documentation.
* Fixed calculation of # 'free' parameters for `aic.mcmc` and `bic.mcmc` criteria when uniquenesses are isotropic.  
  \cr As a result `mixFac_free`, which calculates # 'free' parameters for any _finite_ factor model is exported w/ accompanying documentation.
* Added new plot when `plot.meth=GQ` for OMFA/OMIFA/IMFA/IMIFA methods depicting the trace(s) of the #s of active/non-empty groups.

## Improvements
* Deferred setting of `dimnames` attributes from `mcmc_IMIFA` to `get_IMIFA_results` to reduce memory burden & speed up simulations.
* Jettisoned superfluous duplicate material in the object outputted from `get_IMIFA_results` to reduce size & simplify access.
* Removed IMFA/IMIFA `trunc.G` arg, made `range.G` the max allowable number of active groups & also stored number of active groups.
* Code sped up when G=1 by not simulating labels for OMFA/OMIFA/IMFA/IMIFA & not simulating mixing proportions for OMFA/OMIFA.
* To reduce chance of crash due to exceeding memory capacity, `score.switch` defaults to `FALSE` if number of models being run is large.
* DIC model selection criterion now also available for infinite factor models.
* Invisibly returned from `sim_IMIFA_data`.

## Bug Fixes 
* 2nd label switching move for IMFA/IMIFA sped up & fixed/reweighted to properly ensure uniform sampling of neighbouring pairs of clusters.
* Fixed trace plots for factor scores by extracting indices of stored iterations properly using `Rfast::sort_unique`. 
* Slightly rearranged order in which Gibbs updates take place to ensure means enter simulation of uniquenesses properly.
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA` output.
* Tightened controls for when certain parameters are not stored for posterior inference.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Indices of uncertain observations now returned (`get_IMIFA_results`) & printed (`plot.Results_IMIFA`) even when `zlabels` not supplied.
* Fixed behaviour of progress bar when `verbose=FALSE`.
* Fixed typos and expanded/clarified help documentation & vignette.

# IMIFA v1.1.0 - (__1st release: 2017-02-02__)
