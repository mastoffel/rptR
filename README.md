<!-- README.md is generated from README.Rmd. Please edit that file -->
rptR
====

![Build Status](https://travis-ci.org/mastoffel/rptR.svg?branch=master)

`rptR` provides a collection of functions for calculating point estimates, confidence intervals and significance tests of the repeatability (intra-class correlation coefficient) of measurements, as well as on the variances themselves. The function `rpt` is a the core functions that calls more specialised functions as required. Specialised functions can also be called directly (see `?rpt` for details). All functions return lists of values. The function `?summary.rpt` produces summaries in a detailed format, whereby `?plot.rpt` plots the distributions of bootstrap or permutation test estimates.

-   get the latest development version from github with

``` r
    # install.packages("devtools")
    devtools::install_github("mastoffel/rptR", build_vignettes = TRUE)
    # tutorial
    browseVignettes("rptR")
```
