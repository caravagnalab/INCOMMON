# Getter for class `'INCOMMON'`.

Get model priors used for classification.

## Usage

``` r
priors_k_m(x)
```

## Arguments

- x:

  An obj of class `'INCOMMON'`.

## Value

A dplyr::tibble containing prior distributions.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Get priors used in classification
priors_k_m(MSK_PAAD_output)
#> # A tibble: 26,280 × 6
#> # Groups:   gene, tumor_type [730]
#>    gene  tumor_type     k     m     N     n
#>    <chr> <chr>      <int> <int> <dbl> <dbl>
#>  1 ABL1  HCC            1     1 1026.  55.7
#>  2 ABL1  HCC            2     1 1026. 481. 
#>  3 ABL1  HCC            2     2 1026.  90.6
#>  4 ABL1  HCC            3     1 1026.  87.9
#>  5 ABL1  HCC            3     2 1026.  73.0
#>  6 ABL1  HCC            3     3 1026.  38.9
#>  7 ABL1  HCC            4     1 1026.  25.1
#>  8 ABL1  HCC            4     2 1026.  42.7
#>  9 ABL1  HCC            4     3 1026.  21.7
#> 10 ABL1  HCC            4     4 1026.  16.8
#> # ℹ 26,270 more rows
```
