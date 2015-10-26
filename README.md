<!-- README.md is generated from README.Rmd. Please edit that file -->
rptR
====

![Build Status](https://travis-ci.org/mastoffel/rptR.svg?branch=master)

`rptR` provides a collection of functions for caluculating point estimates, interval estimates and significance tests of the repeatability (intra-class correlation coefficient) of measurements. The function `rpt` is a the core functions that calls more specialised functions as required. Specialised functions can also be called directly (see `?rpt` for details). All functions return lists of values. The function `?summary.rpt` produces summaries in a detailed format.

-   the latest development version from github with

``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("mastoffel/rptR")
    ```
```
