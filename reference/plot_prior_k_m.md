# Plot joint prior distribution of total copy number and multiplicity

Visualises the joint prior distribution of total copy number (\\k\\) and
mutation multiplicity (\\m\\) for each mutation in a sample. The most
likely configuration for each mutation is highlighted.

## Usage

``` r
plot_prior_k_m(priors_k_m, x, k_max)
```

## Arguments

- priors_k_m:

  Prior distribution object for joint \\k,m\\.

- x:

  An INCOMMON object.

- k_max:

  Integer specifying the maximum total copy number.

## Value

A faceted `ggplot2` contour plot showing the joint prior distribution of
copy number and multiplicity for each mutation.
