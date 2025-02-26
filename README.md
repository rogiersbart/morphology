
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The {morphology} R package<img src="man/figures/logo.png" align="right" width="25%"/><br><small><font color="#999">Morphological description of 3D categorical arrays</font></small>

<!-- badges: start -->

[![GitHub R package
version](https://img.shields.io/github/r-package/v/rogiersbart/ti?label=version)](https://github.com/rogiersbart/morphology)
[![CRAN
status](https://www.r-pkg.org/badges/version/morphology)](https://CRAN.R-project.org/package=morphology)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The goal of {morphology} is to …

## Install

You can install the latest version of {morphology} with the following:

``` r
if (!require(pak)) install.packages("pak")
pak::pak("rogiersbart/morphology")
```

## Use

The basic {morphology} workflow looks as follows:

``` r
library(morphology)
description <- my_array |>
  look_at("voxels", of_category = 1) |> 
  look_in(direction = "xyz") |>
  look_for(neighbours = 100, within = 10) |> 
  describe() |> 
  scale_by("inverse neighbourhood") |>
  finalise()
```

See [Get started](articles/morphology.html) for a more comprehensive
overview.

## Note
