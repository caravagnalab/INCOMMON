# Visualise the posterior distribution on (k,m) configurations.

Visualise the posterior distribution on (k,m) configurations.

## Usage

``` r
plot_posterior_k_m(x, k_max = NULL, z_km = NULL)
```

## Arguments

- x:

  An object of class INCOMMON.

- k_max:

  The maximum allowed value of total copyu number.

- z_km:

  The list of posterior distributions over (k,m) configurations.

## Value

An object or a list of objects of class `'ggplot2'`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Plot classification results for a specific sample
x = subset_sample(x = MSK_PAAD_output, sample_list = "P-0000142")
plot_posterior_k_m(x = x)
```
