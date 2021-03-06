
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tvcm

The goal of this package is to provide tools to fit time-varying
coefficient models.

## Installation

``` r
devtools::install_github("inqs909/tvcm")
```

### Windows

You will need to install rtools to allow the `devtools` package to
install packages from a remote repository:
<https://cran.r-project.org/bin/windows/Rtools/>.

Make sure to update R as well.

## Functions

The main functions for this package are the `tvcm`, `binary_vcm`, and
`pois_vcm` functions. They are designed to be similar with other
functions that fit models. Additionally, there are functions to compute
the optimal bandwidth using a cross-validation approach.

## Estimation

The estimates are obtained using either a local least-squares or local
likelihood approach. The varying-coefficients are approximated using a
local-linear model.

## Standard Errors

The standard errors are obtained using a bootstrap method. Due the
correlated nature of longitudinal data, bootstrap standard errors
provide more accurate results. The functions will provide both the
standard errors and point-wise percentiles for a value of a
varying-coefficient functions.

## Plotting

The `plot()` function will generate plots of the time-varying
coefficient model.
