# Get the fraction of alleles with the mutation (FAM) values for a gene and cancer type.

Get the fraction of alleles with the mutation (FAM) values for a gene
and cancer type.

## Usage

``` r
show_FAM(x, tumor_type = NULL, gene = NULL)
```

## Arguments

- x:

  An object of class INCOMMON.

- tumor_type:

  The tumour type identifier.

- gene:

  The gene name..

## Value

A table with the estimated FAM values

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
show_FAM(MSK_PAAD_output, tumor_type = 'PAAD', gene = 'TP53')
#> Joining with `by = join_by(id)`
#> # A tibble: 1,346 × 10
#>    sample   gene  gene_role    NV    DP purity purity_map eta_map      FAM class
#>    <chr>    <chr> <chr>     <int> <int>  <dbl>      <dbl>   <dbl>    <dbl> <chr>
#>  1 P-00001… TP53  TSG          53   671    0.4      0.454    240. 2.50e- 1 Low …
#>  2 P-00015… TP53  TSG         370   703    0.6      0.541    347. 1   e+ 0 High…
#>  3 P-00020… TP53  TSG         295   565    0.5      0.564    180. 7.50e- 1 Bala…
#>  4 P-00022… TP53  TSG          21  1020    0.5      0.162    328. 1.25e- 1 Low …
#>  5 P-00024… TP53  TSG         171   548    0.2      0.333    268. 9.99e- 1 High…
#>  6 P-00024… TP53  TSG         195   601    0.2      0.291    326. 1.00e+ 0 High…
#>  7 P-00024… TP53  TSG          62   534    0.3      0.215    299. 9.99e- 1 High…
#>  8 P-00028… TP53  TSG         170   854    0.5      0.593    160. 2.50e- 1 Bala…
#>  9 P-00029… TP53  TSG         265   447    0.7      0.738    293. 1.00e+ 0 High…
#> 10 P-00029… TP53  TSG         160   722    0.6      0.442    255. 1.51e-10 Low …
#> # ℹ 1,336 more rows
```
