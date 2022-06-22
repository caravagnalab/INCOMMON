
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TAPACLOTH

<!-- badges: start -->
<!-- badges: end -->

The goal of TAPACLOTH is to classify somatic mutations from targeted
panel sequencing as subclonal, clonal heterozygous or clonal LOH.

## Installation

You can install the development version of TAPACLOTH from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/TAPACLOTH")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(TAPACLOTH)
## basic example code
x = dplyr::tibble(sample = "test",
           gene = "test gene",
           nv = 50,
           dp = 100,
           vaf = 0.5,
           purity = 1
)

analyse_sample(
  data = x,
  sample_name = "test",
  alpha_level = 1e-3,
  model = "BetaBinomial",
  rho = 0.01
)
#> 
#> ── test ────────────────────────────────────────────────────────────────────────
#> $fit
#> # A tibble: 1 × 7
#>   sample gene         nv    dp   vaf purity class 
#>   <chr>  <chr>     <dbl> <dbl> <dbl>  <dbl> <chr> 
#> 1 test   test gene    50   100   0.5      1 Clonal
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
