# IMIFA v1.1.1
* __Second release: 2017-02-13__
================================

## New Features

## Improvements
* Deferred setting of various `dimnames` attributes from `mcmc_IMIFA()` to `get_IMIFA_results()` to reduce memory burden & speed up simulations.
* Jettisoned superfluous duplicate material in the object outputted from `get_IMIFA_results()` to reduce size & simplify access.
* Sped up simulation of the mixing proportions for MFA/MIFA/OMFA/OMIFA methods using the normalized `rgamma(G, alpha, 1)` 'trick' from Devroye's book (pg.594) rather than `MCMCpack::rdirichlet`.
* Removed `trunc.G` argument, thereby making `range.G` the maximum allowable number of active groups & also stored the number of active groups for IMFA/IMIFA methods.
* Code sped up when G=1 by avoiding simulating labels for OMFA/OMIFA/IMFA/IMIFA & avoiding simulating mixing proportions for OMFA/OMIFA.
* Sped up simulation of the cluster labels for OMFA/OMIFA/IMFA/IMIFA methods using two part construction & search of the log-cdf.
* Invisibly returned from `sim_IMIFA_data()`.

## Bug Fixes 
* Second label switching move for IMFA/IMIFA methods fixed/reweighted to properly ensure truly uniform sampling of neighbouring pairs of clusters, and also sped up.
* Fixed trace plots for factor scores by extracting indices of stored iterations properly using `Rfast::sort_unique()`. 
* Edited/robustified subsetting of large objects when storing `mcmc_IMIFA()` output.
* Tightened controls for when certain items are not stored.
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.
* Fixed typos and expanded/clarified help documentation & vignette.

# IMIFA v1.1.0 (__1st release: 2017-02-02__)
