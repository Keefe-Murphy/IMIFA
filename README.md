[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/IMIFA)](https://cran.r-project.org/package=IMIFA)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/IMIFA?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/IMIFA?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# IMIFA R Package
## Infinite Mixture of Infinite Factor Analysers
### Written by Keefe Murphy

The IMIFA package provides flexible Bayesian estimation of Infinite Mixtures of Infinite Factor Analysers and related models, for nonparametric model-based clustering of high-dimensional data, introduced by Murphy et al. (2017) \href{https://arxiv.org/abs/1701.07010}{arXiv:1701.07010}. The IMIFA model assumes factor analytic covariance structures within mixture components and simultaneously achieves dimension reduction and clustering without recourse to model selection criteria to choose the number of clusters or cluster-specific latent factors, mostly via efficient Gibbs updates. Model-specific diagnostic tools are also provided, as well as many options for plotting results and conducting posterior inference on parameters of interest.

To install the development version of the package type:

```
# If required install devtools:
# install.packages('devtools')
devtools::install_github('Keefe-Murphy/IMIFA')
```

To install the latest stable official release of the package from CRAN go to R and type:

```
install.packages('IMIFA')
```

In either case, you can then explore the package with:

```
library(IMIFA)
help(mcmc_IMIFA) # Help on the main modelling function
```

To read the vignette guide to using the package, type the following within R:

```
vignette('IMIFA', package = 'IMIFA')
```
