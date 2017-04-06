__Infinite Mixtures of Infinite Factor Analysers__
==================================================
# IMIFA v1.2.0 - (_2nd release: 2017-03-31_)

## New Features
* Learning the Pitman-Yor `discount` and `alpha` parameters via Metropolis-Hastings now implemented.  
  Plotting function's `param` argument gains the option `discount` for posterior inference.
* Sped up simulating cluster labels from unnormalised log probabilities using the Gumbel-Max trick (Yellott, 1977):  
  `gumbel_max` replaces earlier function to sample cluster labels & is now unhidden/exported/documented.
* Added new plot when `plot.meth=GQ` for OM(I)FA/IM(I)FA methods depicting trace(s) of #s of active/non-empty groups.
* Added function `Zsimilarity` to summarise posterior clustering by the sampled labels with minimum  
  squared distance to a sparse similarity matrix constructed by averaging the adjacency matrices.  
  When optionally called inside `get_IMIFA_results`, the similarity matrix can be plotted via `plot.meth="zlabels"`.

## Improvements
* Metropolis-Hastings updates implemented for `alpha` when `discount` is non-zero, rather than usual Gibbs.  
  Mutation rate monitored rather than acceptance rate for Metropolis-Hastings updates of `discount` parameter.
* Fixed calculation of # 'free' parameters for `aic.mcmc` and `bic.mcmc` criteria when uniquenesses are isotropic:    
  `PGMM_dfree`, which calculates # 'free' parameters for _finite_ factor analytic mixture models is exported/documented.  
  This function is also used to add checks on the Dirichlet hyperparameter for OM(I)FA methods.
* DIC model selection criterion now also available for infinite factor models (previously finite only).
* `G_priorDensity` now better reflects discrete nature of the density and plots for non-zero PY discount values.
* Posterior mean loadings heatmaps now also display a colour key legend via new function `heat_legend`.
* Simulated (finite) mixing proportions w/ _Gamma(alpha, 1)_ trick (Devroye 1986, p.594) instead of `MCMCpack:rdirichlet`:  
  `rDirichlet` replaces earlier function to sample mixing proportions & is now unhidden/exported/documented.
* Deferred setting `dimnames` attributes in `mcmc_IMIFA` to `get_IMIFA_results`: lower memory burden/faster simulations.
* Jettisoned superfluous duplicate material in object outputted from `get_IMIFA_results` to reduce size/simplify access.
* IMFA/IMIFA `trunc.G` arg, the max allowable # active groups, defaults to `range.G` & # active groups now stored.
* Code sped up when G=1 by not simulating labels for OM(I)FA/IM(I)FA & not simulating mixing proportions for OM(I)FA.
* Reduced chance of crash by exceeding memory capacity; `score.switch` defaults to `FALSE` if # models ran is large..

## Bug Fixes 
* 2<sup>nd</sup> IM(I)FA label switching move sped up/properly weighted to ensure uniform sampling of neighbouring cluster pairs.
* Offline label switching square assignment correction now permutes properly.
* Fixed factor score trace plots by extracting indices of stored samples using `Rfast::sort_unique` & rotating properly. 
* Fixed adding of `rnorm` columns to scores matrix during adaptation, esp. when widest loadings matrix grows/shrinks.
* Updates of DP/PY `alpha` parameter now correctly depend on current # non-empty rather than active groups.
* Fixed density plots for parameters with bounded support, accounting for spike at zero for `discount`.
* Slightly rearranged order Gibbs updates take place, esp. to ensure means enter simulation of uniquenesses properly.
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA` output.
* Tightened controls for when certain parameters are not stored for posterior inference.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Uniquenesses now stored correctly for all clustering methods.
* Indices of uncertain obs. returned (`get_IMIFA_results`)/printed (`plot.Results_IMIFA`) even when `zlabels` not supplied.
* Fixed behaviour of progress bar when `verbose=FALSE`.
* Fixed typos & expanded/clarified help documentation/vignette.

# IMIFA v1.1.0 - (_1st release: 2017-02-02_)
