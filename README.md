<!-- README.md is generated from README.Rmd. Please edit that file -->
rptR
====

![Build Status](https://travis-ci.org/mastoffel/rptR.svg?branch=master) [![CRAN total downloads](http://cranlogs.r-pkg.org/badges/grand-total/rptR?color=blue)](https://cran.r-project.org/package=rptR)

`rptR` provides a collection of functions for calculating point estimates, confidence intervals and significance tests of the repeatability (intra-class correlation coefficient) of measurements, as well as on the variances themselves. The function `rpt` is a the core functions that calls more specialised functions as required. Specialised functions can also be called directly (see `?rpt` for details). All functions return lists of values. The function `?summary.rpt` produces summaries in a detailed format, whereby `?plot.rpt` plots the distributions of bootstrap or permutation test estimates.

-   get the latest development version from github with

``` r
    # install.packages("devtools")
    # building vignettes might take some time. Set build_vignettes = FALSE for a quick download.
    devtools::install_github("mastoffel/rptR", build_vignettes = TRUE)
    # tutorial
    browseVignettes("rptR")
```

### Reference

Stoffel, M., Nakagawa, S. & Schielzeth, H. (2017) rptR: Repeatability estimation and variance decomposition by generalized linear mixed-effects models. Methods Ecol Evol. Accepted Author Manuscript. <doi:10.1111/2041-210X.12797>
