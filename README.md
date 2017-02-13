[![cran version](http://www.r-pkg.org/badges/version/IMIFA)](https://cran.rstudio.com/web/packages/IMIFA) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/IMIFA?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/IMIFA?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# IMIFA R Package
## Infinite Mixture of Infinite Factor Analysers
### Written by Keefe Murphy

IMIFA is a Bayesian implementation in R of Infinite Mixtures of Infinite Factor Analysers (IMIFA) and related models. The package provides flexible Gibbs sampler functions for fitting IMIFA and related family of models, introduced by Murphy et al. (2017) <https://arxiv.org/abs/1701.07010>, which conducts Bayesian nonparametric model-based clustering with factor analytic covariance structures without recourse to model selection criteria to choose the number of clusters or cluster-specific latent factors. Model-specific diagnostic tools are also provided, as well as many options for plotting results and conducting posterior inference on parameters of interest.

To install the development version of the package type:

```
# If required install devtools:
#install.packages('devtools')
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
