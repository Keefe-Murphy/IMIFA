[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/IMIFA)](https://cran.r-project.org/package=IMIFA)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/IMIFA?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/IMIFA?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# IMIFA R Package
## Infinite Mixture of Infinite Factor Analysers
## (and related models)
### Written by Keefe Murphy

## Description

The IMIFA package provides flexible Bayesian estimation of Infinite Mixtures of Infinite Factor Analysers and related models, for nonparametric model-based clustering of high-dimensional data, introduced by Murphy et al. (2017) <[arXiv:1701.07010](https://arxiv.org/abs/1701.07010)>. The IMIFA model assumes factor analytic covariance structures within mixture components and simultaneously achieves dimension reduction and clustering without recourse to model selection criteria to choose the number of clusters or cluster-specific latent factors, mostly via efficient Gibbs updates. Model-specific diagnostic tools are also provided, as well as many options for plotting results and conducting posterior inference on parameters of interest.

The package also contains two data sets: `olive` and `coffee`.

## Installation

You can install the latest stable official release of the `IMIFA` package from CRAN:

```
install.packages("IMIFA")
```

or the development version from GitHub:

```
# If required install devtools:  
# install.packages('devtools')  
devtools::install_github('Keefe-Murphy/IMIFA')
```

In either case, you can then explore the package with:

```
library(IMIFA)
help(mcmc_IMIFA) # Help on the main modelling function
```

For a more thorough intro, the vignette document is available as follows:

```
vignette("IMIFA", package="IMIFA")
```

However, if the package is installed from GitHub the vignette is not automatically created. It can be accessed when installing from GitHub with the code:

```
devtools::install_github('Keefe-Murphy/IMIFA', build_vignettes = TRUE)
```

Alternatively, the vignette is available on the package's CRAN page.

### References
Murphy, K., Gormley, I. C. and Viroli, C. (2017) Infinite Mixtures of Infinite Factor Analysers: Nonparametric Model-Based Clustering via Latent Gaussian Models. _To appear_. <[arXiv:1701.07010](https://arxiv.org/abs/1701.07010)>
