__Infinite Mixtures of Infinite Factor Analysers__
==================================================  

### New Features
* `get_IMIFA_results` gains the arg. `adapt`, allowing automatic, post-hoc truncation of  
redundant factors in cases where `mcmc_IMIFA` was called with `adapt=FALSE` only  
(provided relevant parameters were stored) for the `"IFA"`, `"MIFA"`, `"OMIFA"`, & `"IMIFA"` methods.
* `mgpControl` gains the `active.crit` arg. to allow the criterion for determining the number of  
active factors of Schiavon & Canale (2020) to be used (for `"IFA"` models only); this is available  
for both `mcmc_IMIFA` with `adapt=TRUE` directly and the aforementioned post-hoc adaptation using  
`get_IMIFA_results`, though the old criterion of Bhattacharya & Dunson (2011), `active.crit="BD"`,  
remains the default for all infinite factor models in each case. See documentation for details. 

### Improvements, Bug Fixes, & Miscellaneous Edits
* Models with `range.Q=0` now work again for the `"FA"`, `"MFA"`, `"OMFA"`, & `"IMFA"` methods.  
_(with thanks to @abcoxyzide for flagging the issue)_
* The `mgpControl` arg. `delta0g` now also governs `beta.d1`, `beta.d2`,  
  & `sigma.hyper`, in addition to `alpha.d1` & `alpha.d2`.
* Improved documentation & checks related to the `Ledermann` bound + added various references to,  
  & warnings for, the stricter `Q <= floor((P - 1)/2)` restriction throughout the package.
* Minor fixes for `"MFA"`/`"MIFA"` models when `mu0g`, `psi0g`, &/or `delta0g`  
  is `TRUE`, particularly when `length(range.G) > 1`.
* Minor speed-up to `G_priorDensity` when `discount` is non-zero.
* Extensive edits to avoid overheads introduced in `matrixStats` v1.0.0 + related minor speed-ups.
* Now using newer `CITATION` commands.

## IMIFA v2.1.10 - (_17<sup>th</sup> release [patch update]: 2022-12-19_)
### Improvements, Bug Fixes, & Miscellaneous Edits
* Minor fix to adaptive Gibbs sampler for clusters with zero factors.
* Minor speed-ups to `Procrustes` when `translate` &/or `dilate` are `TRUE`.
* Minor speed-ups to internal Cholesky decompositions in rare failure cases.
* Cosmetic changes to `param="alpha"` & `param="discount"` plots.
* Minor speed-ups, documentation clarifications, & examples edits.

## IMIFA v2.1.9 - (_16<sup>th</sup> release [patch update]: 2022-08-12_)
### Improvements, Bug Fixes, & Miscellaneous Edits
* `storeControl` gains the `update.mu` arg. to optionally circumvent updates of the  
  mean parameters to speed-up special case of (I)FA models with centered data.
* Minor fix to `param="means"` plots for uncentered (I)FA results.
* Minor fixes to account for (rare) uncentered data problems in means updates also.
* Minor fix to rare cases where Cholesky decompositions fail.
* Minor tidy-ups related to `exchange=TRUE` from previous update.
* Minor documentation edits for CRAN compliance.

## IMIFA v2.1.8 - (_15<sup>th</sup> release [patch update]: 2021-12-19_)
### Improvements, Bug Fixes, & Miscellaneous Edits
* Slight modifications to adaptive Gibbs sampling for all infinite factor methods:  
  * Adaptation now occurs before related parameter updates: cleaner, slightly faster code.
  * Fixes __rare__ bug re: storing factor scores when there are empty components.
  * Associated new default behaviour(s) for `mgpControl` arg. `start.AGS`.
* `bnpControl` gains the args. `thresh` & `exchange`, implementing the exchangeable/thresholded  
  slice sampler of Fall & Barat (2014): these are experimental features (both args. default to `FALSE`).
* Fixed posterior predictive checking bug in `get_IMIFA_results` for  
  models for univariate data where _some_ components have zero factors.
* Minor IM(I)FA speed-ups to updates of DP concentration parameter alpha.
* `sapply` replaced with `vapply`, with other negligible speed-ups.

## IMIFA v2.1.7 - (_14<sup>th</sup> release [patch update]: 2021-10-07_)
### Improvements, Bug Fixes, & Miscellaneous Edits
* `mgpControl` gains the arg. `truncated` (defaults to `FALSE`):  
  * Allows version of MGP prior with gamma distributions left-truncated at 1  
  to be used for the column shrinkage parameters beyond the first column:  
  this has some more desirable shrinkage properties, at the expense of longer run times.
  * Related function `rltrgamma` to draw from left-truncated gamma distributions provided.
  * `exp_ltrgamma` for computing expectations of such distributions also provided, and used  
  internally by `MGP_check` when its own new `truncated` arg. is `TRUE` (default: `FALSE`).
  * See `?ltrgamma` for further details on `rltrgamma` and `exp_ltrgamma`.
* Minor speed-ups to simulation of component mean parameters.
* Major speed-up to simulation of factor loadings parameters (especially when Q=1).
* Major speed-up to simulation of factor scores when Q=1.
* Minor speed-ups to simulation of means and loadings from priors for empty components.
* Minor fixes to adaptive Gibbs sampler for MIFA models when `delta0g=TRUE`.
* Improved checks on `range.G` & `range.Q` in `mcmc_IMIFA`.
* Minor fixes to returned attributes from `MGP_check` when invoked in vectorised fashion.
* Minor vignette styling edits and documentation clarifications.

## IMIFA v2.1.6 - (_13<sup>th</sup> release [patch update]: 2021-05-24_)
### Bug Fixes & Miscellaneous Edits
* Fixed breaking bugs associated with IM(I)FA slice samplers introduced in previous update.
* `G_calibrate` function exported to augment existing `G_expected` & `G_variance` functions.
* `G_variance` now computed more accurately and efficiently for the `alpha=0` case.
* Major speed-up to `G_expected` for the `alpha=0` case.
* Minor speed-ups to simulation of local/column/cluster shrinkage parameters + scores & loadings.
* Minor speed-up to `G_priorDensity` for non-zero `discount`.
* Minor speed-up to `psi_hyper`.

## IMIFA v2.1.5 - (_12<sup>th</sup> release [patch update]: 2020-12-29_)
### Bug Fixes & Miscellaneous Edits
* Fixed `mcmc_IMIFA` initialisation issues for univariate data sets.
* `error.metrics=TRUE` now works for univariate data sets in `get_IMIFA_results`.
* Slight speed-ups to slice samplers for IM(I)FA methods.
* Prettier axis labels for first two plots produced by `plot.meth="zlabels"`.
* Added ORCID iDs to DESCRIPTION.
* Minor CRAN compliance edits to the vignette.

## IMIFA v2.1.4 - (_11<sup>th</sup> release [patch update]: 2020-11-18_)
### Bug Fixes & Miscellaneous Edits
* Stricter checking of permissible `alpha` values in the special case of negative `discount`.
* Fixes when `plot.meth="density"` and `param="alpha"` for _fixed_ positive `discount`.
* Cosmetic changes to `G_priorDensity` plots.
* Ensured `matrixStats (>= 0.53.1)` and `mclust (>= 5.4)` in `Imports:`.
* Ensured `gmp (>= 0.5-4)` in `Suggests:`.
* Package startup message now checks if newer version of package is available from CRAN.
* Minor documentation & examples edits (esp. `mat2cols`).
* Replaced dependency on `viridis` with `viridisLite`.
* Updated citation info after publication in _Bayesian Analysis_.
* Updated maintainer e-mail address.

## IMIFA v2.1.3 - (_10<sup>th</sup> release [patch update]: 2020-05-12_)
### Bug Fixes & Miscellaneous Edits
* Maintenance release for compatibility with R 4.0.0 - minor edits.
* Improved handling of suggested packages `Rmpfr` & `gmp` in `G_expected`, `G_variance`, & `G_priorDensity`.
* `summary.Results_IMIFA` gains the printing-related argument `MAP=TRUE`.
* Edited printed details when `plot.meth="zlabels"` with unsupplied `zlabels`.
* Minor fixes for fixed negative `discount` (an experimental feature).
* Minor speed-up to `Procrustes` when `dilate=TRUE` (never used internally).
* Minor efficiency gains in slice samplers for IM(I)FA methods.
* Documentation, vignette, examples, and references improvements.

## IMIFA v2.1.2 - (_9<sup>th</sup> release [patch update]: 2020-03-30_)
### Bug Fixes
* Fixes and speed-ups to MGP updates and adaptive Gibbs sampler for IMIFA/OMIFA/MIFA models:  
    * Fixes and speed-ups to MGP parameter updates when _some_ clusters have zero factors.
    * Additional speed-ups to simulation of column-shrinkage parameters when _some_ clusters are empty.
    * Properly accounted for the cluster-shrinkage parameters when the number of factors increases.
    * Minor bug fixes for padding scores when the maximum number of factors increases.
* Variable-specific communalities (`x$Error$Var.Exps`) now returned by `get_IMIFA_results` in addition  
  to proportion of explained variance per cluster (`x$Error$Clust.Exps`; previously `x$Error$Var.Exps`).
* `G_expected` & `G_variance` gain the arg. `MPFR` to control use of suggested packages.
* Minor speed-up to `rDirichlet` for the symmetric uniform case.
* Ensured compatibility with latest version of `Rfast` package (w/ minor speed-ups).
* Removed `methods` package from `Suggests:`.
* Spell-checking of documentation and fixes to `donttest` examples.

## IMIFA v2.1.1 - (_8<sup>th</sup> release [patch update]: 2019-12-11_)
### Improvements
* `discount` can now be fixed at a negative value when `learn.d=FALSE`,  
  provided `alpha` is supplied as a positive integer multiple of `abs(discount)` and `learn.alpha=TRUE`.
* Other types of `norm` (beyond Frobenius) can now be specified, by passing the arg. `type`,  
  via the `...` construct, for calculating the PPRE within `get_IMIFA_results`.
* The breaks used to construct the bins for the PPRE calculation can now also be specified,  
  by passing the `dbreaks` arg. through the `...` construct. This is an experimental feature; caution is advised.
* The settings `discount<0` & `alpha=0` now accommodated by `G_expected`, `G_variance`, & `G_priorDensity`:  
  `G_expected` no longer requires the `Rmpfr` or `gmp` libraries for non-zero `discount` unless `alpha=0`. 
* `mgpControl` gains the arg. `forceQg` (defaults to `FALSE`, i.e. retains old behaviour - see documentation for details).
* `G_priorDensity` gains `type` arg. and now works again in non-vectorised form.
* Minor speed-up to `Procrustes` function and hence the identifiability corrections within `get_IMIFA_results`.
* Minor speed-ups to `post_conf_mat` function and `"parallel.coords"` plots.
* Updated citation info after online publication in _Bayesian Analysis_.

### Bug fixes
* Fixes to `sim_IMIFA_data` to allow empty clusters and related fix for `nonempty` arg. to `get_IMIFA_results`.
* Fixed bug when initial `alpha` value is `0` when `learn.alpha=TRUE`.
* Minor fix for handling optional args. to `mixfaControl` and `plot.Results_IMIFA` functions.
* Admissible `rho` values in `bnpControl` corrected to [0,1) from (0,1].
* Fixed bug related to Procrustes rotation of the factor scores for (I)/FA models in `get_IMIFA_results`.
* Fixed handling of colour palettes in `plot.Results_IMIFA` & `G_priorDensity`.
* Documentation and warning message fixes.
* Anti-aliasing of vignette images.

## IMIFA v2.1.0 - (_7<sup>th</sup> release [minor update]: 2019-02-04_)
### New Features
* `mgpControl` gains the arguments `cluster.shrink` and `sigma.hyper`:
    * `cluster.shrink` governs invocation of cluster shrinkage MGP hyperparameter for MIFA/OMIFA/IMIFA methods.
    * `sigma.hyper` controls the gamma hyperprior on this parameter. The posterior mean is reported, where applicable.
    * Full conditionals for loadings and local/column shrinkage MGP hyperparameters edited accordingly.
* Allowed the Dirichlet concentration parameter `alpha` to be learned via MH steps for the OM(I)FA models.  
    * Also allowed diminishing adaptation to tune the log-normal proposal to achieve a target acceptance rate.
    * Thus `bnpControl` args. `learn.alpha`, `alpha.hyper`, `zeta`, & `tune.zeta` become relevant for OM(I)FA models.
* New posterior predictive model checking approach added to `get_IMIFA_results` (with associated plots):  
  Posterior Predictive Reconstruction Error (PPRE) compares bin counts of the original data with corresponding  
  counts for replicate draws from the posterior predictive distribution using a standardised Frobenius norm.
* Added new function `scores_MAP` to decompose factor scores summaries  
  from `get_IMIFA_resuls` into submatrices corresponding to the MAP partition.
* Added new wrapper function `sim_IMIFA_model` to call `sim_IMIFA_data` using  
  the estimated parameters from fitted `Results_IMIFA` objects.
* New `get_IMIFA_results` arg. `vari.rot` allows loadings templates to be varimax rotated,  
  prior to Procrustes rotation, for more interpretable solutions (defaults to `FALSE`).
* New `plot.Results_IMIFA` argument `common` governing `plot.meth="means"` plots (details in documentation).
  
### Improvements
* New hyperparameter/argument defaults:
    * `sigma.mu` defaults to `1` s.t. the hypercovariance is the identity for the prior on the means;  
      old behaviour (using the diagonal of the sample covariance matrix) recoverable by specifying `sigma.mu=NULL`.
    * `prec.mu` defaults to `0.01` s.t. the prior on the cluster means is flat by default.
    * `learn.d` defaults to `TRUE` s.t. a PYP prior is assumed for IM(I)FA models by default.
    * `alpha.hyper` now has a larger hyper-rate by default, to better encourage clustering.
    * `alpha.d1` & `alpha.d2` now set to `2.1`/`3.1` rather than `2`/`6` to discourage exponentially fast shrinkage.
    * `z.init` now defaults to `"hc"`: model-based agglomerative hierarchical clustering.
* Overhauled `psi_hyper` (details in documentation) for:  
    * `N <= P` data where the sample covariance matrix is not invertible.
    * `type="isotropic"` uniquenesses.
* Args. `scores` & `loadings` can now be supplied to `sim_IMIFA_data` directly;  
  new arg. `non.zero` controls the # effective factors (per column & cluster) when `loadings` are instead simulated.
* Sped-up 2<sup>nd</sup> label-switching move for IM(I)FA models (accounting for empty clusters).
* Args. for `hc` can now be passed when `init.z="mclust"` also  
  (previously only `"hc"`), thus controlling how `Mclust` is itself initialised.
* Allowed `criterion` to be passed via `...` in `mixfaControl` to choose between  
  `mclustBIC`/`mclustICL` to determine optimum model to initialise with when  
  `z.init="mclust"` & also sped-up `mclust` initialisation in the process.
* Added `stop.AGS` arg. to `mgpControl`: renamed `adapt.at` to `start.AGS` for consistency.
* Added `start.zeta` & `stop.zeta` options to `tune.zeta` argument in `bnpControl`.
* Allowed user-supplied `breaks` in the plotting functions `mat2cols` & `heat_legend`.
* Initial cluster sizes are now shown in order to alert users to potentially bad starting values.
* Added utility function `pareto_scale()`.

### Bug Fixes
* Fixed factor _scores_ & error metrics issues in `get_IMIFA_results` for clustering methods:  
    * Fixed storage of scores for infinite factor methods - now corresponds to samples where the  
      largest cluster-specific number of factors is `>=` the max of the modal estimates of the same  
      (previously samples where __any__ cluster has `>=` the corresponding modal estimate were used):  
      thus, valid samples for computing error metrics also fixed and Procrustes rotation also sped-up.  
    * Other Procrustes rotation fixes to account for label-switching.  
    * Other Procrustes rotation fixes specific to the IMFA/OMFA methods.
* `range.G` and `trunc.G` defaults fixed, especially for small sample size settings. 
* Slight label-switching fixes when `zlabels` are supplied to `get_IMIFA_results`;  
  posterior confusion matrix, cluster sizes vector, and the sampled labels themselves effected.
* Prevented unnecessary Procrustes rotation for single-factor components, thus fixing some bugs.
* Fixed initialisation of uniquenesses to account for all four settings of `uni.type`.
* Allowed conditioning on iterations with all components populated for M(I)FA models in `get_IMIFA_results`.
* Accounted for 1-component IM(I)FA/OM(I)FA models in `get_IMIFA_results`.
* Fixed handling of empty components when simulating cluster labels from priors in `mcmc_IMIFA` & `sim_IMIFA_data`.
* Ensured no. of factors `Q` cannot exceed no. of observations in the corresponding cluster in `sim_IMIFA_data`.
* Slight speed-up to updating MGP hyperparameters in the presence of empty MIFA/OMIFA/IMIFA components.
* Slight speed-up to sampling cluster labels with slice indicators for IM(I)FA models.
* Explicitly allowed Pitman-Yor special case where `alpha=0` for IM(I)FA models;  
  added related controls on spike-and-slab prior for `discount` when fixing `alpha<=0`.
* Allowed full range of `hc` model types for initialisation purposes via `...` in `mixfaControl`.
* Clarified `dimnames` of `get_IMIFA_results` output in `x$Loadings` & `x$Scores`.
* Fixed storage switches & iteration indices to better account for `burnin=0`.
* Fixed plotting of exact zeros in posterior confusion matrix.
* Fixed plotting posterior mean loadings heatmap when one or more clusters have zero factors.
* Fixed plotting scores for (I)FA models due to bug in previous update, esp. with `zlabels` supplied.
* Fixed `show_IMIFA_digit` to better account for missing pixels &/or the data having been centered/scaled.
* Fixed simulation of `psi` when not supplied to `sim_IMIFA_data` to IG rather than GA.
* Fixed bug preventing `Q` to be supplied to `get_IMIFA_results` for infinite factor methods.
* Fixed y-axis labelling in uncertainty type plots when `plot.meth="zlabels"`.
* Small fixes to function `show_digit`.
* Better handling of tied model-selection criteria in `get_IMIFA_results`.
* `Procrustes` now works when `X` has fewer columns than `Xstar`.
* Minor cosmetic change for overplotting `scores` & `loadings` in `trace` & `density` plots.
* Edited `Ledermann` and related warnings to account for case of isotropic uniquenesses.
* Tidied indentation/line-breaks for `cat`/`message`/`warning` calls for printing clarity.
* Corrected `IMIFA-package` help file (formerly just `IMIFA`).
* Edited `CITATION` file and authorship.

## IMIFA v2.0.0 - (_6<sup>th</sup> release [major update]: 2018-05-01_)
### Major Changes
* Simplified `mcmc_IMIFA` by consolidating arguments using new helper functions (with defaults):  
    * Args. common to all factor-analytic mixture methods & MCMC settings supplied via `mixfaControl`.
    * MGP & AGS args. supplied via `mgpControl` for infinite factor models.  
    * Pitman-Yor/Dirichlet Process args. supplied via `bnpControl` for infinite mixture models.
    * Storage switch args. supplied via `storeControl`.
    * New functions also inherit the old documentation for their arguments.

### New Features
* Posterior predictive checking overhauled: now MSE, RMSE etc. between empirical & estimated covariance  
  matrices are computed for every retained iteration so uncertainty in these estimates can be quantified:  
    * Can be switched on/off via the `error.metrics` argument to `get_IMIFA_results`.  
    * Can be visualised by supplying `plot.meth="errors"` to `plot.Results_IMIFA`.  
    * For methods which achieve clustering, the 'overall' covariance matrix  
      is now properly computed from the cluster-specific covariance matrices.  
    * Same metrics also evaluated at posterior mean parameter estimates & for final sample where possible.
* `mixfaControl` gains the arg. `prec.mu` to control the degree of flatness of the prior for the means.
* Posterior confusion matrix now returned (`get_IMIFA_results`) & visualisable (`plot.Results_IMIFA`,  
  when `plot.meth="zlabels"`), via new function `post_conf_mat`, to further assess clustering uncertainty.
* Added new type of clustering uncertainty profile plot in `plot.Results_IMIFA` when `plot.meth="zlabels"`.
* For convenience, `get_IMIFA_results` now also returns the last valid samples for parameters of interest,  
  after conditioning on the modal G & Q values and accounting for label switching and Procrustes rotation.
* `plot.Results_IMIFA` gains new arg. `show.last` that replaces any instance of showing the posterior mean  
  with the last valid sample instead (i.e. when `plot.meth="means"` or `plot.meth="parallel.coords")`.
* Added ability to constrain mixing proportions across clusters using `equal.pro` argument for M(I)FA models:  
  Modified `PGMM_dfree` accordingly and forced non-storage of mixing proportions when `equal.pro` is `TRUE`. 
* All methods now work for univariate data also (with apt. edits to plots & uniqueness defaults etc.).  
  `sim_IMIFA_data` also extended to work for univariate data, as well as sped-up.

### Improvements
* Retired args. `nu` & `nuplus1` to `mgpControl`, replaced by ability to specify more general gamma prior,  
  via new `phi.hyper` arg. specifying shape _and_ rate - `MGP_check` has also been modified accordingly.
* `Zsimilarity` sped-up via the `comp.psm` & `cltoSim` functions s.t. when # observations < 1000.
* Matrix of posterior cluster membership probabilities now returned by `get_IMIFA_results`.
* Modified AGS to better account for when the number of group-specific latent factors shrinks to zero.
* `psi.alpha` no longer needs to be strictly greater than 1, unless the default `psi.beta` is invoked;  
  thus flatter inverse gamma priors can now be specified for the uniquenesses via `mixfaControl`.
* Added "`hc`" option to `z.init` to initialise allocations via hierarchical clustering (using `mclust::hc`).
* Allowed optional args. for functions used to initialise allocations via `...` in `mixfaControl`.
* Added `mu` argument to `sim_IMIFA_data` to allow supplying true mean parameter values directly.
* Standard deviation of `aicm`/`bicm` model selection criteria now computed and returned.
* Speed-ups due to new `Rfast` utility functions: `colTabulate` & `matrnorm`.
* Speed-ups due to utility functions from `matrixStats`, on which `IMIFA` already depends.
* Slight improvements when `adapt=FALSE` for infinite factor models with fixed high truncation level.
* Misclassified observations now highlighted in 1<sup>st</sup> type of uncertainty plot in `Plot.Results_IMIFA`,  
  when `plot.meth="zlabels"` and the true `zlabels` are supplied.
* `mixfaControl` gains arg. `drop0sd` to control removal of zero-variance features (defaults to `TRUE`).
* `heat_legend` gains `cex.lab` argument to control magnification of legend text.
* `mat2cols` gains the `transparency` argument.
*  Edited `PGMM_dfree` to include the 4 extra models from the EPGMM family.

### Bug Fixes
* Supplying `zlabels` to `get_IMIFA_results` will now match the cluster labels and parameters to  
  the true labels even if there is a mismatch between the number of clusters in both.
* Similarly, supplying `zlabels` to `plot.Results_IMIFA` when `plot.meth="zlabels"` no longer does  
  any matching when printing performance metrics to the screen - previously this caused confusion  
  as associated parameters were not also permuted as they are within `get_IMIFA_results`: now  
  `plot(get_IMIFA_results(sim), plot.meth="zlabels", zlabels=z)` gives different results from  
  `plot(get_IMIFA_results(sim, zlabels=z), plot.meth="zlabels")` as only the latter will permute.
* Accounted for errors in covariance matrix when deriving default `sigma.mu` & `psi.beta` values.
* Accounted for missing empirical covariance entries within `get_IMIFA_results`.
* Fixed model selection in `get_IMIFA_results` for IMFA/OMFA models when `range.Q` is a range.
* Fixed calculation of `aicm`, `bicm` and `dic` criteria: all results remain the same.
* Fixed support of Ga(a, b) prior on `alpha` when `discount` is being learned.
* Fixed bug preventing `uni.prior="isotropic"` when `uni.type` is `(un)constrained`.
* Fixed treatment of exact zeros when plotting average clustering similarity matrix.
* Fixed tiny bug when neither centering nor scaling (of any kind) are applied to data within `mcmc_IMIFA`.
* Fixed plotting of posterior mean scores when one or more clusters are empty.
* Fixed bug with default plotting palette for data sets with >1024 variables.
* Fixed bug with label switching permutations in `get_IMIFA_results` when there are empty clusters.
* Fixed `print` and `summary` functions for objects of class `IMIFA` and `Results_IMIFA`.
* Fixed calculating posterior mean `zeta` when adaptively targeting `alpha`'s optimal MH acceptance rate.
* Allowed `alpha` be tiny for (O)M(I)FA models (provided `z.init != "priors"` for overfitted models).
* Normalised mixing proportions in `get_IMIFA_results` when conditioning on `G` for IM(I)FA/OM(I)FA models.
* New controls/warnings for excessively small Gamma hyperparemeters for uniqueness/local shrinkage priors.
* Clarified recommendation in `MGP_check` that `alpha.d2` be moderately large relative to `alpha.d1`.
* Ensured `sigma.mu` hyperparameter arg. is always coerced to diagonal entries of a covariance matrix.
* Transparency default in `plot.Results_IMIFA` now depends on device's support of semi-transparency.
* Replaced certain instances of `is.list(x)` with `inherits(x, "list")` for stricter checking.
* Added `check.margin=FALSE` to calls to `sweep`.
* `Ledermann`, `MGP_check`, & `PGMM_dfree` are now properly vectorised.

### Miscellaneous Edits
* Added `USPSdigits` data set (training and test),  
  with associated utility functions `show_digit` and `show_IMIFA_digit`.
* Optimised compression of `olive`, `coffee` and vignette data and used `LazyData: true`.
* Added `call.=FALSE` to `stop()` messages and `immediate.=TRUE` to certain `warning()` calls.
* Removed dependency on`adrop`, `e1071`, `graphics`, `grDevices`, `plotrix`, `stats` & `utils` libraries.
* Reduced dependency on `Rfast` w/ own version of `standardise`.
* Added utility function `IMIFA_news` for accessing this `NEWS` file.
* Added `CITATION` file.
* Extensively improved package documentation: 
    * Added `Collate:` field to `DESCRIPTION` file.
    * Added line-breaks to `usage` sections of multi-argument functions.
    * Consolidated help files for `G_expected` & `G_variance`.

## IMIFA v1.3.1 - (_5<sup>th</sup> release [patch update]: 2017-07-07_)
### Bug Fixes
* Fixed bug preventing M(I)FA models from being treated as (I)FA models when `range.G` contains 1.
* Fixed bug preventing `get_IMIFA_results` from working properly when true labels are NOT supplied.

## IMIFA v1.3.0 - (_4<sup>th</sup> release [minor update]: 2017-06-22_)
### New Features
* Added options `"constrained"` & `"single"` to `mcmc_IMIFA`'s `uni.type` argument:  
  as well as being either diagonal or isotropic (UUU / UUC), uniquenesses can now further be  
  constrained across clusters (UCU / UCC), with appropriate warnings, defaults, checks,  
  initialisations, computation of model choice penalties, and plotting behaviour in all 4 cases.
* `mcmc_IMIFA` gains the `tune.zeta` argument, a list of `heat`, `lambda` & `target` parameters, to invoke  
  diminishing adaptation for tuning the uniform proposal to achieve a target acceptance rate when `alpha`  
  is learned via  Metropolis-Hastings when the Pitman-Yor Process prior is employed for the IM(I)FA models.
  
### Improvements
* (I)FA models sped up by considering uniquenesses under 1-cluster models as `"constrained"` or `"single"`,  
  rather than previously `"unconstrained"` or `"isotropic"`, utilising pre-computation and empty assignment.
* Previously hidden functions improved, exported and documented with examples:  
  `is.cols`, `Ledermann`, `Procrustes` & `shift_GA`.
* `is.posi_def` gains `make` argument, merging it with previously hidden function `.make_posdef`:  
  Thus the 'nearest' positive-(semi)definite matrix and the usual check can be returned in a single call.
* Sped-up sampling IM(I)FA labels, esp. when 'active' G falls to 1, or the _dependent_ slice-sampler is used:  
  `log.like` arg. removed from `gumbel_max`; function stands alone, now only stored log-likelihoods computed.
* `psi` argument added to `sim_IMIFA_data` to allow supplying true uniqueness parameter values directly.
  
### Bug Fixes
* Used `bw="SJ"` everywhere `density` is invoked for plotting (`bw="nrd0"` is invoked if this fails).
* Fixed initialisation of uniquenesses for `isotropic` (I)FA models.
* Fixed parallel coordinates plot axes and labels for all `isotropic` uniquenesses plots.
* Fixed adaptation for MIFA/OMIFA/IMIFA models when all clusters simultaneously have zero factors.
* Fixed storage bug in IM(I)FA models when `learn.d` is `TRUE` but `learn.alpha` is `FALSE`.
* Fixed density plot for `discount` when mutation rate is too low (i.e. too many zeros).
* Fixed simulation of loadings matrices for empty MIFA/OMIFA/IMIFA clusters using `byrow=TRUE`:  
  Loop to simulate loadings matrices now generally faster also for all models.
* Fixed silly error re: way in which (I)FA models are treated as 1-cluster models to ensure they run:  
  Related bug fixed for OM(I)FA/IM(I)FA models when starting number of clusters is actually supplied.

## IMIFA v1.2.1 - (_3<sup>rd</sup> release [patch update]: 2017-05-29_)
### Improvements
* Posterior mean scores can now also be plotted in the form of a heat map (previously loadings only).  
  `load.meth` argument replaced by logical `heat.map` in `plot.Results_IMIFA`.
* `mat2cols` gains `compare` argument to yield common palettes/breaks for heat maps of multiple matrices:  
  Associated `plot_cols` function also fixed, and now unhidden.
* Removed certain dependencies with faster personal code: e.g. Procrustes rotation now quicker:  
  `IMIFA` no longer depends on the `corpcor`, `gclus`, `MASS`, `matrixcalc`, or `MCMCpack` libraries.
  
### Bug Fixes 
* Used `par()$bg` (i.e. default `"white"`) for plotting zero-valued entries of similarity matrix.
* Range of data for labelling in `heat_legend` calculated correctly.
* `mcmc_IMIFA`'s `verbose` argument now governs printing of `message` & `cat` calls, but not `stop` or `warning`.
* Fixed storage and plotting of loadings, particularly when some but not all clusters have zero factors.
* Added `NEWS.md` to build.

## IMIFA v1.2.0 - (_2<sup>nd</sup> release [minor update]: 2017-05-09_)
### New Features
* Learning the Pitman-Yor `discount` & `alpha` parameters via Metropolis-Hastings now implemented.  
    * Spike-and-slab prior specified for `discount`: size of spike controlled by arg. `kappa`.
    * Plotting function's `param` argument gains the option `discount` for posterior inference.
* Sped up simulating cluster labels from unnormalised log probabilities using the Gumbel-Max trick (Yellott, 1977):  
  `gumbel_max` replaces earlier function to sample cluster labels and is now unhidden/exported/documented.
* Added new plot when `plot.meth=GQ` for OM(I)FA/IM(I)FA models depicting trace of #s of active/non-empty clusters.
* Added function `Zsimilarity` to summarise posterior clustering by the sampled labels with minimum  
  squared distance to a sparse similarity matrix constructed by averaging the adjacency matrices:  
  when optionally called inside `get_IMIFA_results`, the similarity matrix can be plotted via `plot.meth="zlabels"`.

### Improvements
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
* Restored the IMFA/IMIFA arg. `trunc.G`, the max allowable # active clusters, and # active clusters now stored.
* Code sped up when `active` G=1 by not simulating labels for IM(I)FA models.
* Reduced chance of crash by exceeding memory capacity; `score.switch` defaults to `FALSE` if # models ran is large.

### Bug Fixes 
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

## IMIFA v1.1.0 - (_1<sup>st</sup> release: 2017-02-02_)
