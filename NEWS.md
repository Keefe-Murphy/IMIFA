__Infinite Mixtures of Infinite Factor Analysers__
==================================================  

# IMIFA v1.3.2 - (_6<sup>th</sup> release [patch update]: 2017-07-26_)
## Improvements
* Added ability to constrain mixing proportions across clusters using `equal.pro` argument for M(I)FA models:  
  modified PGMM_dfree accordingly and forced non-storage of mixing proportions when `equal.pro` is TRUE.  
* All methods now work for univariate data also (with apt. edits to plots & uniqueness defaults etc.).
* `Zsimilarity` sped-up via `mcclust::comp.psm` s.t. `z.avgsim=TRUE` now by default in `get_IMIFA_results`.
* Added "`hc`" option to `z.init` to initialise allocations via hierarchical clustering (using `mclust::hc`).
* Simplified `mcmc_IMIFA` by consolidating storage switch arguments using new function `storeControl`.
* Standard deviation of AICM/BICM model selection criteria now computed and returned.
* `heat_legend` gains `cex.lab` argument to control magnification of legend text.

## Bug Fixes
* Fixed bug preventing `uni.prior="isotropic"` when `uni.type` is `(un)constrained`.
* Fixed model selection in `get_IMIFA_results` for IMFA/OMFA models when `range.Q` is a range.
* Fixed treatment of exact zeros when plotting average clustering similarity matrix.
* Fixed tiny bug when neither centering nor scaling (of any kind) are applied to data within `mcmc_IMIFA`.
* Fixed plotting of posterior mean scores when one or more clusters are empty.
* Fixed storage switches to account for `burnin=0`.
* Fixed bug with default plotting palette for data sets with >1024 variables.
* Transparency default in `plot.Results_IMIFA` now depends on device's support of semi-transparency.

# IMIFA v1.3.1 - (_5<sup>th</sup> release [patch update]: 2017-07-07_)
## Bug Fixes
* Fixed bug preventing M(I)FA models from being treated as (I)FA models when `range.G` contains 1.
* Fixed bug preventing `get_IMIFA_results` from working properly when true labels are NOT supplied.

# IMIFA v1.3.0 - (_4<sup>th</sup> release [minor update]: 2017-06-22_)
## New Features
* Added options `"constrained"` & `"single"` to `mcmc_IMIFA`'s `uni.type` argument:  
  as well as being either diagonal or isotropic (UUU / UUC), uniquenesses can now further be  
  constrained across clusters (UCU / UCC), with appropriate warnings, defaults, checks,  
  initialisations, computation of model choice penalties, and plotting behaviour in all 4 cases.
* `mcmc_IMIFA` gains the `tune.zeta` argument, a list of `heat`, `lambda` & `target` parameters, to invoke  
  diminishing adaptation for tuning the uniform proposal to achieve a target acceptance rate when `alpha`  
  is learned via  Metropolis-Hastings when the Pitman-Yor Process prior is employed for the IM(I)FA models.
  
## Improvements
* (I)FA models sped up by considering uniquenesses under 1-cluster models as `"constrained"` or `"single"`,  
  rather than previously `"unconstrained"` or `"isotropic"`, utilising pre-computation and empty assignment.
* Previously hidden functions improved, exported and documented with examples:  
  `is.cols`, `Ledermann`, `Procrustes` & `shift_GA`.
* `is.posi_def` gains `make` argument, merging it with previously hidden function `.make_posdef`:  
  Thus the 'nearest' positive-(semi)definite matrix and the usual check can be returned in a single call.
* Sped-up sampling IM(I)FA labels, esp. when 'active' G falls to 1, or the _dependent_ slice-sampler is used:  
  `log.like` arg. removed from `gumbel_max`; function stands alone, now only stored log-likelihoods computed.
* `psi` argument added to `sim_IMIFA_data` to allow supplying true uniqueness parameter values directly.
  
## Bug Fixes
* Used `bw="SJ"` everywhere `density` is invoked for plotting (`bw="nrd0"` is invoked if this fails).
* Fixed initialisation of uniquenesses for `isotropic` (I)FA models.
* Fixed parallel coordinates plot axes and labels for all `isotropic` uniquenesses plots.
* Fixed adaptation for MIFA/OMIFA/IMIFA models when all clusters simultaneously have zero factors.
* Fixed storage bug in IM(I)FA models when `learn.d` is `TRUE` but `learn.alpha` is `FALSE`.
* Fixed density plot for `discount` when mutation rate is too low (i.e. too many zeros).
* Fixed simulation of loadings matrices for empty MIFA/OMIFA/IMIFA clusters using `byrow=TRUE`:  
  loop to simulate loadings matrices now generally faster also for all models.
* Fixed silly error re: way in which (I)FA models are treated as 1-cluster models to ensure they run:  
  Related bug fixed for OM(I)FA/IM(I)FA models when starting number of clusters is actually supplied.

# IMIFA v1.2.1 - (_3<sup>rd</sup> release [patch update]: 2017-05-29_)
## Improvements
* Posterior mean scores can now also be plotted in the form of a heat map (previously loadings only).  
  `load.meth` argument replaced by logical `heat.map` in `plot.Results_IMIFA`.
* `mat2cols` gains `compare` argument to yield common palettes/breaks for heat maps of multiple matrices:  
  Associated `plot_cols` function also fixed, and now unhidden.
* Removed certain dependencies with faster personal code: e.g. Procrustes rotation now quicker:  
  `IMIFA` no longer depends on the `corpcor`, `gclus`, `MASS`, `matrixcalc`, or `MCMCpack` libraries.
  
## Bug Fixes 
* Used `par()$bg` (i.e. default `"white"`) for plotting zero-valued entries of similarity matrix.
* Range of data for labelling in `heat_legend` calculated correctly.
* `mcmc_IMIFA`'s `verbose` argument now governs printing of `message` & `cat` calls, but not `stop` or `warning`.
* Fixed storage and plotting of loadings, particularly when some but not all clusters have zero factors.
* Added `NEWS.md` to build.

# IMIFA v1.2.0 - (_2<sup>nd</sup> release [minor update]: 2017-05-09_)

## New Features
* Learning the Pitman-Yor `discount` & `alpha` parameters via Metropolis-Hastings now implemented.  
  Plotting function's `param` argument gains the option `discount` for posterior inference.
* Sped up simulating cluster labels from unnormalised log probabilities using the Gumbel-Max trick (Yellott, 1977):  
  `gumbel_max` replaces earlier function to sample cluster labels and is now unhidden/exported/documented.
* Added new plot when `plot.meth=GQ` for OM(I)FA/IM(I)FA models depicting trace of #s of active/non-empty clusters.
* Added function `Zsimilarity` to summarise posterior clustering by the sampled labels with minimum  
  squared distance to a sparse similarity matrix constructed by averaging the adjacency matrices:  
  when optionally called inside `get_IMIFA_results`, the similarity matrix can be plotted via `plot.meth="zlabels"`.

## Improvements
* Metropolis-Hastings updates implemented for `alpha` when `discount` is non-zero, rather than usual Gibbs.  
  Mutation rate monitored rather than acceptance rate for Metropolis-Hastings updates of `discount` parameter.
* Fixed calculation of # '_free_' parameters for `aic.mcmc` & `bic.mcmc` criteria when uniquenesses are isotropic:    
  `PGMM_dfree`, which calculates # 'free' parameters for _finite_ factor analytic mixture models is exported/documented.  
  This function is also used to add checks on the Dirichlet hyperparameter for OM(I)FA methods.
* DIC model selection criterion now also available for infinite factor models (previously finite only).
* `G_priorDensity` now better reflects discrete nature of the density, and plots for non-zero PY discount values.
* Posterior mean loadings heatmaps now also display a colour key legend via new function `heat_legend`.
* Avoided redundant simulation of stick-breaking/mixing proportions under both types of IM(I)FA slice sampler.
* Simulated (finite) mixing proportions w/ _Gamma(alpha, 1)_ trick (Devroye 1986, p.594) instead of `MCMCpack:rdirichlet`:  
  `rDirichlet` replaces earlier function to sample mixing proportions and is now unhidden/exported/documented.
* Deferred setting `dimnames` attributes in `mcmc_IMIFA` to `get_IMIFA_results`: lower memory burden/faster simulations.
* Jettisoned superfluous duplicate material in object outputted from `get_IMIFA_results` to reduce size/simplify access.
* IMFA/IMIFA `trunc.G` arg, the max allowable # active clusters, defaults to `range.G` and # active clusters now stored.
* Code sped up when `active` G=1 by not simulating labels for IM(I)FA models.
* Reduced chance of crash by exceeding memory capacity; `score.switch` defaults to `FALSE` if # models ran is large.

## Bug Fixes 
* 2<sup>nd</sup> IM(I)FA label switching move sped up/properly weighted to ensure uniform sampling of neighbouring cluster pairs.
* Offline label switching square assignment correction now permutes properly.
* Fixed factor score trace plots by extracting indices of stored samples using `Rfast::sort_unique` and rotating properly. 
* Fixed adding of `rnorm` columns to scores matrix during adaptation, esp. when widest loadings matrix grows/shrinks.
* Fixed initialisation (and upper limit) of number of clusters for OM(I)FA/IM(I)FA, esp. when `N < P`.
* Updates of DP/PY `alpha` parameter now correctly depend on current # non-empty rather than active clusters.
* Fixed density plots for parameters with bounded support, accounting for spike at zero for `discount`.
* Slightly rearranged order Gibbs updates take place, esp. to ensure means enter simulation of uniquenesses properly.
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA` output.
* Tightened controls for when certain parameters are not stored for posterior inference.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Geometric rather than arithmetic mean used to derive single rate hyperparameter for PPCA's isotropic uniquenesses.
* Uniquenesses now stored correctly for all clustering methods.
* Indices of uncertain obs. returned (`get_IMIFA_results`)/printed (`plot.Results_IMIFA`) even when `zlabels` not supplied.
* Fixed behaviour of progress bar when `verbose=FALSE`.
* Fixed typos and expanded/clarified help documentation/vignette.

# IMIFA v1.1.0 - (_1<sup>st</sup> release: 2017-02-02_)
