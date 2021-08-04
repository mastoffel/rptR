<!-- README.md is generated from README.Rmd. Please edit that file -->

# rptR

![Build Status](https://travis-ci.org/mastoffel/rptR.svg?branch=master)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/rptR?color=blue)](https://cran.r-project.org/package=rptR)
[![Codecov test
coverage](https://codecov.io/gh/mastoffel/rptR/branch/master/graph/badge.svg)](https://codecov.io/gh/mastoffel/rptR?branch=master)

### Goal

The goal of `rptR` is to provide point estimates, confidence intervals
and significance tests of the **repeatability** (intra-class correlation
coefficient) of measurements based on generalised linear mixed models
(GLMMs). The function `?summary.rpt` produces summaries in a detailed
format, whereby `?plot.rpt` plots the distributions of bootstrap or
permutation test estimates.

### Citation

When using \`rptR\`\`, please cite our
[paper](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12797):

Stoffel, M. A., Nakagawa, S., & Schielzeth, H. (2017). rptR:
Repeatability estimation and variance decomposition by generalized
linear mixed-effects models. Methods in Ecology and Evolution, 8(11),
1639-1644.
[![](https://img.shields.io/badge/doi-https://doi.org/10.1111/2041--210X.12797-green.svg)](https://doi.org/https://doi.org/10.1111/2041-210X.12797)

### Installation

You can install the stable version of \`rptR\`\` from CRAN with:

``` r
install.packages("rptR")
```

Or the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("mastoffel/rptR", build_vignettes = TRUE, dependencies = TRUE) 
# manual
browseVignettes("rptR")
```

### Example

``` r
library(rptR)
data(BeetlesBody)

rpts <- rpt(BodyL ~ Treatment + Sex + (1 | Container) + (1 | Population), 
            grname = c("Container", "Population"), data = BeetlesBody, 
            datatype = "Gaussian", nboot = 100, npermut = 100)
#> Bootstrap Progress:
#> Permutation Progress for Container :
#> Permutation Progress for Population :
```

``` r
summary(rpts)
#> 
#> Repeatability estimation using the lmm method
#> 
#> Call = rpt(formula = BodyL ~ Treatment + Sex + (1 | Container) + (1 | Population), grname = c("Container", "Population"), data = BeetlesBody, datatype = "Gaussian", nboot = 100, npermut = 100)
#> 
#> Data: 960 observations
#> ----------------------------------------
#> 
#> Container (120 groups)
#> 
#> Repeatability estimation overview: 
#>       R     SE   2.5%  97.5% P_permut  LRT_P
#>  0.0834  0.025 0.0452  0.132     0.01      0
#> 
#> Bootstrapping and Permutation test: 
#>             N    Mean   Median   2.5%  97.5%
#> boot      100 0.08772 0.084018 0.0452 0.1321
#> permut    100 0.00463 0.000405 0.0000 0.0238
#> 
#> Likelihood ratio test: 
#> logLik full model = -1528.553
#> logLik red. model = -1555.264
#> D  = 53.4, df = 1, P = 1.34e-13
#> 
#> ----------------------------------------
#> 
#> 
#> Population (12 groups)
#> 
#> Repeatability estimation overview: 
#>       R     SE   2.5%  97.5% P_permut  LRT_P
#>   0.491 0.0978  0.294  0.652     0.04      0
#> 
#> Bootstrapping and Permutation test: 
#>             N   Mean Median   2.5%  97.5%
#> boot      100  0.477  0.472  0.294  0.652
#> permut    100  0.452  0.455  0.419  0.492
#> 
#> Likelihood ratio test: 
#> logLik full model = -1528.553
#> logLik red. model = -1595.399
#> D  = 134, df = 1, P = 3.19e-31
#> 
#> ----------------------------------------
```

#### Plotting

``` r
plot(rpts, grname="Container", type="boot", cex.main=0.8, color = "#5E81AC")
#> Warning in plot.window(xlim, ylim, "", ...): "color" is not a graphical
#> parameter
#> Warning in title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...): "color"
#> is not a graphical parameter
#> Warning in axis(1, ...): "color" is not a graphical parameter
#> Warning in axis(2, ...): "color" is not a graphical parameter
plot(rpts, grname="Population", type="boot", cex.main=0.8, color = "#5E81AC")
#> Warning in plot.window(xlim, ylim, "", ...): "color" is not a graphical
#> parameter
#> Warning in title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...): "color"
#> is not a graphical parameter
#> Warning in axis(1, ...): "color" is not a graphical parameter
#> Warning in axis(2, ...): "color" is not a graphical parameter
```

![](README-gaussian-1.png)![](README-gaussian-2.png)
