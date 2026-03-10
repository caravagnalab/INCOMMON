# Visualise the prior distribution on sample purity

Visualise the prior distribution on sample purity

## Usage

``` r
plot_purity_prior(x, sample, purity_error = 0.05)
```

## Arguments

- x:

  An object of class INCOMMON.

- sample:

  An identifier of the tumour sample.

- purity_error:

  The variance of the Beta prior distribution on purity.

## Value

An object or a list of objects of class `'ggplot2'`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Plot classification results for a specific sample
plot_purity_prior(x = MSK_PAAD_output, sample = "P-0000142", purity_error = 0.05)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the INCOMMON package.
#>   Please report the issue at <https://github.com/caravagnalab/INCOMMON/issues>.
```
