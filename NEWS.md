# IMIFA v1.1.1 - (__2nd release: 2017-02-15__)
============================================

## New Features
* Sped up simulating cluster labels from unnormalised log probabilities using the Gumbel-Max trick (Yellott, 1977):  
  `gumbel_max` replaces earlier function to sample cluster labels & is now unhidden/exported/documented.
* Simulated (finite) mixing proportions w/ _Gamma(alpha, 1)_ trick (Devroye 1986, p.594) instead of `MCMCpack:rdirichlet`:  
  `rDirichlet` replaces earlier function to sample mixing proportions & is now unhidden/exported/documented.
* Fixed calculation of # 'free' parameters for `aic.mcmc` and `bic.mcmc` criteria when uniquenesses are isotropic:    
  `PGMM_dfree`, which calculates # 'free' parameters for _finite_ factor analytic mixture models is exported/documented.  
  This function is also used to add checks on the Dirichlet hyperparameter for OM(I)FA methods.
* Added new plot when `plot.meth=GQ` for OM(I)FA/IM(I)FA methods depicting trace(s) of #s of active/non-empty groups.

## Improvements
* Metropolis-Hastings updates implemented for `alpha` when `discount` is non-zero, rathern than Gibbs.
* Deferred setting `dimnames` attributes in `mcmc_IMIFA` to `get_IMIFA_results`: lower memory burden/faster simulations.
* Jettisoned superfluous duplicate material in object outputted from `get_IMIFA_results` to reduce size/simplify access.
* Removed IMFA/IMIFA `trunc.G` arg, made `range.G` the max allowable # active groups & also stored # active groups.
* Code sped up when G=1 by not simulating labels for OM(I)FA/IM(I)FA & not simulating mixing proportions for OM(I)FA.
* Reduced chance of crash by exceeding memory capacity; `score.switch` defaults to `FALSE` if # models ran is large.
* DIC model selection criterion now also available for infinite factor models (previously finite only).
* `G_priorDensity` now better reflects discrete nature of the density and plots for non-zero PY discount values.
* Invisibly returned from `sim_IMIFA_data`.

## Bug Fixes 
* 2<sup>nd</sup> IM(I)FA label switching move sped up/properly weighted to ensure uniform sampling of neighbouring cluster pairs.
* Fixed trace plots for factor scores by extracting indices of stored iterations properly using `Rfast::sort_unique`. 
* Fixed way in which `rnorm` columns are added to scores matrix during adaptation when 'widest' loadings matrix grows.
* Slightly rearranged order Gibbs updates take place, esp. to ensure means enter simulation of uniquenesses properly.
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA` output.
* Tightened controls for when certain parameters are not stored for posterior inference.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Uniquenesses now stored correctly for all clustering methods.
* Indices of uncertain obs. returned (`get_IMIFA_results`)/printed (`plot.Results_IMIFA`) even when `zlabels` not supplied.
* Fixed behaviour of progress bar when `verbose=FALSE`.
* Fixed typos & expanded/clarified help documentation/vignette.

# IMIFA v1.1.0 - (__1st release: 2017-02-02__)
