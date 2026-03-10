# Gene roles from COSMIC Cancer Gene Census

This table contains gene roles in cancer, as obtained from the COSMIC
Cancer Gene Census v98. Data were curated such that each gene in list is
assigned either TSG or oncogene.

## Usage

``` r
data(cancer_gene_census)
```

## Format

A data frame with 733 rows and 2 columns:

- gene:

  Name of the gene (Hugo Symbol)

- gene_role:

  Tumour Suppressor Gene (TSG) or oncogene

## Source

COSMIC Cancer Gene Census: https://cancer.sanger.ac.uk/census

## Examples

``` r
data(cancer_gene_census)
cancer_gene_census
#> # A tibble: 733 × 2
#>    gene   gene_role
#>    <chr>  <chr>    
#>  1 A1CF   oncogene 
#>  2 ABI1   TSG      
#>  3 ABL1   oncogene 
#>  4 ABL2   oncogene 
#>  5 ACKR3  oncogene 
#>  6 ACSL3  oncogene 
#>  7 ACSL6  oncogene 
#>  8 ACVR1  oncogene 
#>  9 ACVR1B TSG      
#> 10 ACVR2A TSG      
#> # ℹ 723 more rows
```
