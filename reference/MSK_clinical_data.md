# Genomic data of the MSK-MetTropism cohort

This table contains sample names and corresponding clinical data from
the MSK-MetTropism cohort.

## Usage

``` r
data(MSK_clinical_data)
```

## Format

A tibble with 25659 rows and 15 columns:

## Source

MSK-MET at cBioPortal
https://www.cbioportal.org/study/summary?id=msk_met_2021

## Examples

``` r
data(MSK_clinical_data)
MSK_clinical_data
#> # A tibble: 25,368 × 17
#> # Rowwise: 
#>    sample    tumor_type purity OS_MONTHS OS_STATUS SAMPLE_TYPE MET_COUNT
#>    <chr>     <chr>       <dbl>     <dbl>     <dbl> <chr>           <dbl>
#>  1 P-0000004 BRCA          0.5      3.78         1 Primary             2
#>  2 P-0000015 BRCA          0.4     13.9          1 Metastasis          8
#>  3 P-0000024 UCEC          0.4     35.1          1 Metastasis          8
#>  4 P-0000025 UCEC          0.3     46            1 Metastasis         13
#>  5 P-0000026 UCEC          0.1     80.6          0 Metastasis         11
#>  6 P-0000034 BLCA          0.4      0.79         1 Primary             4
#>  7 P-0000037 HCC           0.9     80.9          0 Metastasis          1
#>  8 P-0000039 PLEMESO       0.4      5.62         1 Primary             5
#>  9 P-0000041 BRCA          0.3     13.6          1 Primary             9
#> 10 P-0000042 PLEMESO       0.4     56.9          1 Primary             0
#> # ℹ 25,358 more rows
#> # ℹ 10 more variables: METASTATIC_SITE <chr>, MET_SITE_COUNT <dbl>,
#> #   PRIMARY_SITE <chr>, SUBTYPE_ABBREVIATION <chr>, GENE_PANEL <chr>,
#> #   SEX <chr>, TMB_NONSYNONYMOUS <dbl>, FGA <dbl>, AGE_AT_SEQUENCING <dbl>,
#> #   RACE <chr>
```
