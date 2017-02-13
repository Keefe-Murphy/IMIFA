# IMIFA v1.1.1
* __Second release: 2017-02-13__
================================

## New Features

## Improvements
* Invisibly returned from `sim_IMIFA_data()`.
* Jettisoned superfluous duplicate material in the object outputted from `get_IMIFA_results()` to reduce size and simplify access.
* Deferred setting of various `dimnames` attributes from `mcmc_IMIFA()` to `get_IMIFA_results()` to reduce memory burden and speed up simulations.
* Sped-up simulation of the mixing proportions for MFA/MIFA/OMFA/OMIFA methods using the normalized `rgamma(G, alpha, 1)` 'trick' from Devroye's book (pg.594) rather than `MCMCpack::rdirichlet`.
* Code sped up when G=1 by avoiding simulating labels for OMFA/OMIFA/IMFA/IMIFA and avoiding simulating mixing proportions for OMFA/OMIFA.
* Removed `trunc.G` argument and stored the number of active groups for IMFA/IMIFA methods.

## Bug Fixes 
* Fixed typos and expanded/clarified help documentation.
* Tightened controls for when certain items are not stored.
* Second label switching move for IMFA/IMIFA methods fixed/reweighted to properly ensure truly uniform sampling and also sped-up.
* Edited subsetting or large objects when storing `mcmc_IMIFA()` output.
* Fixed trace plots for factor scores by extracting indices of stored iterations properly using `Rfast::sort_unique()`. 
* Edited Ledermann upper bound `stop(...)` for finite factor models to `warning(...)`.

# IMIFA v1.1.0
* __First release: 2017-02-02__
================================
