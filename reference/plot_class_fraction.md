# Visualize frequency distribution of INCOMMON classes.

Visualize frequency distribution of INCOMMON classes.

## Usage

``` r
plot_class_fraction(x, tumor_type = NULL, gene = NULL, ...)
```

## Arguments

- x:

  A list of objects of class `'INCOMMON'` containing the classification
  results, as produced by using function `classify`.

- tumor_type:

  Tumor type for tumor-specific prior ('PANCA' for pan-cancer).

- gene:

  Gene for gene-specific prior.

- ...:

  Default S3 method parameter.

## Value

An object or a list of class `'ggplot2'`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Plot class fraction for a specific gene and tumour type
MSK_PAAD_output = mutant_dosage_classification(MSK_PAAD_output)
#> Joining with `by = join_by(id)`
plot_class_fraction(x = MSK_PAAD_output, tumor_type = 'PAAD', gene = 'KRAS')
```
