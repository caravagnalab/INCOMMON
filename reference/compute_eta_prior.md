# Evaluate a prior distribution over the rate of reads per chromosom copy from the data.

Evaluate a prior distribution over the rate of reads per chromosom copy
from the data.

## Usage

``` r
compute_eta_prior(x, priors_k_m = priors_k_m)
```

## Arguments

- x:

  An object of class INCOMMON.

- priors_k_m:

  Pre-computed priors for gene muation total copy number and
  multiplicity.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
data(priors_pcawg_hmf)
compute_eta_prior(x = MSK_PAAD_output, priors_k_m = priors_pcawg_hmf)
#> # A tibble: 2 × 6
#>   tumor_type mean_eta var_eta     N alpha_eta beta_eta
#>   <chr>         <dbl>   <dbl> <int>     <dbl>    <dbl>
#> 1 PAAD           336.  17039.  1659      6.62   0.0197
#> 2 PANCA          336.  17039.  1659      6.62   0.0197
```
