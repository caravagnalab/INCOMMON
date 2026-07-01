# Generate a diagnostic report of an INCOMMON model fit

Combines, for a single sample, the posterior predictive checks for the
Poisson (depth) and Binomial (variant reads) sub-models, the eta and
purity posterior predictive checks, and the prior/posterior
distributions of copy number and multiplicity into one composite figure.

## Usage

``` r
plot_report(x)
```

## Arguments

- x:

  A classified object of class INCOMMON (see
  [`classify`](caravagnalab.github.io/INCOMMON/reference/classify.md)).

## Value

A composite `ggplot2`/`patchwork` object with the full diagnostic report
for the sample.

## Examples

``` r
if (FALSE) { # \dontrun{
# plot_report requires a classified INCOMMON object (see ?classify)
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
plot_report(x)
} # }
```
