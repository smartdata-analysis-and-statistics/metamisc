# metamisc: an R package for meta-analysis of prognosis studies

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/metamisc)](https://cran.r-project.org/package=metamisc)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=metamisc)
<!-- badges: end -->

This is the official repository of the R package **metamisc**, which was developed to facilitate meta-analysis of diagnosis and prognosis research studies. The R package includes functions to  summarize multiple estimates of prediction model discrimination and calibration performance (Debray et al., 2019) <doi:10.1177/0962280218785504>. It also includes functions to evaluate funnel plot asymmetry (Debray et al., 2018) <doi:10.1002/jrsm.1266>. Finally, the package provides functions for developing multivariable prediction models from datasets with clustering (de Jong et al., 2021) <doi:10.1002/sim.8981>

## Installation

The `metamisc` package can be installed from CRAN as follows:

``` r
install.packages("metamisc")
```

The latest version can be installed from GitHub as follows:

``` r
install.packages("devtools")
devtools::install_github(repo = "smartdata-analysis-and-statistics/metamisc")
```

## JASP
A visual interface to the software has been implemented by JASP <https://jasp-stats.org/>
