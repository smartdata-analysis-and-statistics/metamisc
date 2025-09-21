
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metamisc <img src="man/figures/logo.png" align="right" width="150" alt="" />

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/metamisc)](https://cran.r-project.org/package=metamisc)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/metamisc)](https://cran.r-project.org/package=metamisc)
[![R-CMD-check](https://github.com/smartdata-analysis-and-statistics/metamisc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/smartdata-analysis-and-statistics/metamisc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is the official repository of the R package **metamisc**, which was
developed to facilitate meta-analysis of diagnosis and prognosis
research studies. The package includes functions for the following
tasks:

- To develop and validate multivariable prediction models from datasets
  with clustering (de Jong et al., 2021)
- To summarize multiple estimates of prediction model discrimination and
  calibration performance (Debray et al., 2019)
- To evaluate funnel plot asymmetry (Debray et al., 2018)

## Installation

The `metamisc` package can be installed from CRAN as follows:

``` r
install.packages("metamisc")
```

You can install the development version of metamisc from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("smartdata-analysis-and-statistics/metamisc")
```

## JASP

A visual interface to the software has been implemented by JASP
<https://jasp-stats.org/>

## Funding

The development of this R package has been funded by the following
organisations:

- The Netherlands Organisation for Health Research and Development
  (grant 91617050).
- The European Union’s Horizon 2020 research and innovation programme
  under ReCoDID grant agreement No 825746.
- Smart Data Analysis and Statistics B.V., a limited liability
  corporation registered at the Netherlands Chamber of Commerce under
  number 863595327.

## References

de Jong VMT, Moons KGM, Eijkemans MJC, Riley RD, Debray TPA. Developing
more generalizable prediction models from pooled studies and large
clustered data sets. *Stat Med*. 2021 May 5;40(15):3533–59.

Debray TPA, Moons KGM, Riley RD. Detecting small-study effects and
funnel plot asymmetry in meta-analysis of survival data: a comparison of
new and existing tests. *Res Syn Meth*. 2018;9(1):41–50.

Debray TPA, Damen JAAG, Riley R, Snell KIE, Reitsma JB, Hooft L, et
al. A framework for meta-analysis of prediction model studies with
binary and time-to-event outcomes. *Stat Methods Med Res*. 2019
Sep;28(9):2768–86.
