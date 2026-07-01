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

## Examples

``` r
if (FALSE) { # \dontrun{
# plot_prior_k_m requires a classified INCOMMON object (see ?classify)
data(MSK_genomic_data)
data(MSK_clinical_data)
data(priors_pcawg_hmf)
data(priors_eta)
sample = 'P-0002081'
x = init(
  genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,],
  clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,]
)
x = classify(
  x = x, priors_k_m = priors_pcawg_hmf, priors_eta = priors_eta,
  num_cores = 1, iter_warmup = 10, iter_sampling = 10, num_chains = 1
)
plot_prior_k_m(priors_k_m = x$priors_k_m, x = x, k_max = 8)
} # }
```
