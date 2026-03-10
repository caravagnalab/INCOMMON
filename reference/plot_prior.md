# Visualize prior distribution for a gene (tumor-specific or pancancer).

Visualize prior distribution for a gene (tumor-specific or pancancer).

## Usage

``` r
plot_prior(x, gene, tumor_type)
```

## Arguments

- x:

  A prior distribution in the format required for `INCOMMON`, such as
  [`INCOMMON::priors_pcawg_hmf`](caravagnalab.github.io/INCOMMON/reference/priors_pcawg_hmf.md).

- gene:

  Gene for gene-specific prior.

- tumor_type:

  Tumor type for tumor-specific prior ('PANCA' for pan-cancer).

## Value

An object or a list of objects of class `'ggplot2'`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Plot classification results for a specific sample
plot_prior(x = MSK_PAAD_output, gene = 'TP53', tumor_type = 'PAAD')
```
