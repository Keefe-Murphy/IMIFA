# IMIFA v1.1.1 - (__2nd release: 2017-02-15__)
============================================

## New Features
* Sped up simulating cluster labels for OMFA/OMIFA/IMFA/IMIFA methods using 2-part construction & search of necessary parts of the log-cdf. \cr As a result `sim_z_log` replaces earlier functions to simulate cluster labels & is now unhidden & exported with accompanying documentation.
* Added new plot when `plot.meth=GQ` for OMFA/OMIFA/IMFA/IMIFA methods depicting the trace(s) of the #s of active/non-empty groups.

## Improvements
* Deferred setting of `dimnames` attributes from `mcmc_IMIFA` to `get_IMIFA_results` to reduce memory burden & speed up simulations.
* Jettisoned superfluous duplicate material in the object outputted from `get_IMIFA_results` to reduce size & simplify access.
* Simulated finite/overfitted mixing proportions using `rgamma(G, alpha, 1)` trick (Devroye, p.594) instead of `MCMCpack:rdirichlet`.
* Removed IMFA/IMIFA `trunc.G` arg, made `range.G` the max allowable number of active groups & also stored number of active groups.
* Code sped up when G=1 by not simulating labels for OMFA/OMIFA/IMFA/IMIFA & not simulating mixing proportions for OMFA/OMIFA.
* Indices of uncertain observations now returned (`get_IMIFA_results`) & printed (`plot.Results_IMIFA`) even when `zlabels` not supplied.
* Invisibly returned from `sim_IMIFA_data`.

## Bug Fixes 
* 2nd label switching move for IMFA/IMIFA sped up & fixed/reweighted to properly ensure uniform sampling of neighbouring pairs of clusters.
* Fixed trace plots for factor scores by extracting indices of stored iterations properly using `Rfast::sort_unique`. 
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA` output.
* Tightened controls for when certain paramters are not stored for posterior inference.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Slightly rearranged order in which Gibbs updates take place.
* Fixed behaviour of progress bar when `verbose=FALSE`.
* Fixed typos and expanded/clarified help documentation & vignette.

# IMIFA v1.1.0 - (__1st release: 2017-02-02__)
