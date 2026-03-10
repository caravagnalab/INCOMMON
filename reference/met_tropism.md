# Logistic regression of metastatic tropism based on INCOMMON classes.

Logistic regression of metastatic tropism based on INCOMMON classes.

## Usage

``` r
met_tropism(x, gene, tumor_type, metastatic_site)
```

## Arguments

- x:

  An object of class `'INCOMMON'` containing the classification results
  as produced by function `classify`.

- gene:

  The gene on which patient's stratification is based.

- tumor_type:

  The tumor type of patients to stratify.

- metastatic_site:

  The target organ of metastatic diffusion.

## Value

An object of class `'INCOMMON'` containing an additional object
`survival`. Logistic regression of metastatic tropism based on INCOMMON
classes.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Estimate the metastatic propensity associated with mutant TP53 with vs
# without CNA in BRCA, with respect to the Liver.
MSK_PAAD_output = met_tropism(x = MSK_PAAD_output, tumor_type = 'BRCA', gene = 'TP53', metastatic_site = 'Liver')
#> Joining with `by = join_by(id)`
```
