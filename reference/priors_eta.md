# Priors over the rate of reads per chromosome copy from MSK-MET

Prior distribution of copy number and mutation multiplicity from PCAWG.

## Usage

``` r
data(priors_eta)
```

## Format

A data frame with 12 rows and 6 columns:

- tumor_type:

  Tumor type.

- mean_eta:

  Mean of the distribution

- var_eta:

  Variance of the distribution

- N:

  Total number of samples

- alpha_eta:

  Shape parameter alpha of the distribution

- beta_eta:

  Shape parameter beta of the distribution

## Source

Validated copy number calls from PCAWG:
https://doi.org/10.5281/zenodo.6410935

## Examples

``` r
data(priors_eta)
priors_eta
#> # A tibble: 12 × 6
#>    tumor_type mean_eta var_eta     N alpha_eta beta_eta
#>    <chr>         <dbl>   <dbl> <int>     <dbl>    <dbl>
#>  1 BRCA           320.  31113.  2271      3.28  0.0103 
#>  2 ESCA           357.  28769.   285      4.44  0.0124 
#>  3 GIST           340.  61998.   299      1.87  0.00549
#>  4 HCC            304.  16885.   171      5.48  0.0180 
#>  5 LUAD           353.  34458.  3406      3.61  0.0102 
#>  6 MEL            304.  28999.  1042      3.19  0.0105 
#>  7 OV             316.  19945.   936      5.01  0.0159 
#>  8 PAAD           334.  15101.  1698      7.38  0.0221 
#>  9 PRAD           283.  21340.   223      3.75  0.0133 
#> 10 SCLC           372.  45882.   255      3.01  0.00811
#> 11 STAD           298.  17672.    85      5.04  0.0169 
#> 12 PANCA          332.  29352. 10732      3.75  0.0113 
```
