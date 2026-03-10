# Getter for class `'INCOMMON'`.

Get model parameters of the performed classification tests.

## Usage

``` r
parameters(x)
```

## Arguments

- x:

  An obj of class `'INCOMMON'`.

## Value

A dplyr::tibble containing parameters for all the models used in the
classification.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Get classification parameters
parameters(MSK_PAAD_output)
#> # A tibble: 1 × 11
#>   k_max purity_error num_cores iter_warmup iter_sampling num_chains results_dir 
#>   <dbl>        <dbl>     <dbl>       <dbl>         <dbl>      <dbl> <chr>       
#> 1     8         0.05         4        1000          2000          4 ~/INCOMMON_…
#> # ℹ 4 more variables: generate_report_plot <lgl>, reports_dir <chr>,
#> #   stan_fit_dump <lgl>, stan_fit_dir <chr>
```
