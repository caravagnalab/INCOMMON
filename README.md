
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TAPACLOTH <a href="caravagnalab.github.io/TAPACLOTH"><img src="man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/TAPACLOTH/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/caravagnalab/TAPACLOTH/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

TAPACLOTH is a tool that classifies somatic mutation calls from
targeted-panel sequencing as either subclonal, clonal heterozygous or
clonal LOH, thorugh a statistical hypothesis test based on a Binomial or
Beta-Binomial modeling of the sequencing read count process.
Additionally, TAPACLOTH can be used to infer sample purity. Since the
two processes are independent, one can first estimate sample purity and
then use it as input for the classification task.

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/TAPACLOTH/-yellow.svg)](https://caravagnalab.github.io/TAPACLOTH)

## Installation

You can install the development version of TAPACLOTH from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/TAPACLOTH")
```

## Example

This is a basic example of running the classification:

``` r
library(TAPACLOTH)
## basic example code
x = dplyr::tibble(sample = "test",
           gene = "test gene",
           nv = 50,
           dp = 100,
           VAF = 0.5,
           purity = 1
)

run_classifier(
  data = x,
  alpha_level = 1e-3,
  model = "BetaBinomial",
  rho = 0.01
)
#> 
#> ── test gene ───────────────────────────────────────────────────────────────────
#> $fit
#> # A tibble: 1 × 7
#>   sample gene         nv    dp   VAF purity class_terzile
#>   <chr>  <chr>     <dbl> <dbl> <dbl>  <dbl> <chr>        
#> 1 test   test gene    50   100   0.5      1 Clonal LOH   
#> 
#> $model
#> [1] "BetaBinomial"
#> 
#> $rho
#> [1] 0.01
#> 
#> $alpha_level
#> [1] 0.001
```
