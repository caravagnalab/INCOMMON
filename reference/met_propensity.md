# Logistic regression of metastatic propensity based on INCOMMON classes.

Logistic regression of metastatic propensity based on INCOMMON classes.

## Usage

``` r
met_propensity(x, gene, tumor_type)
```

## Arguments

- x:

  An object of class `'INCOMMON'` containing the classification results
  as produced by function `classify`.

- gene:

  The gene on which patient's stratification is based.

- tumor_type:

  The tumor type of patients to stratify.

## Value

An object of class `'INCOMMON'` containing an additional object
`survival`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Estimate the metastatic propensity associated with mutant TP53 with vs without CNA in BRCA.
MSK_PAAD_output = met_propensity(x = MSK_PAAD_output, tumor_type = 'BRCA', gene = 'TP53')
#> Joining with `by = join_by(id)`
```
