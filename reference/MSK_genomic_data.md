# Genomic data of the MSK-MetTropism cohort

This table contains sample names and corresponding read counts data from
the MSK-MetTropism cohort.

## Usage

``` r
data(MSK_genomic_data)
```

## Format

A tibble with 224939 rows and 10 columns:

## Source

MSK-MET at cBioPortal
https://www.cbioportal.org/study/summary?id=msk_met_2021

## Examples

``` r
data(MSK_genomic_data)
MSK_genomic_data
#> # A tibble: 224,939 × 10
#>    sample    chr        from        to ref   alt      DP    NV   VAF gene   
#>    <chr>     <chr>     <dbl>     <dbl> <chr> <chr> <int> <int> <dbl> <chr>  
#>  1 P-0028912 chr17   7577121   7577121 G     A       837   133 0.159 TP53   
#>  2 P-0028912 chr6  111983080 111983081 -     A       698   141 0.202 FYN    
#>  3 P-0028912 chrX   53246994  53246994 G     A       832    85 0.102 KDM5C  
#>  4 P-0003698 chr17   7576852   7576852 C     A       437   109 0.249 TP53   
#>  5 P-0003698 chr3   49933259  49933259 C     A       591    86 0.146 MST1R  
#>  6 P-0003698 chr5  149435631 149435631 C     T       360    36 0.1   CSF1R  
#>  7 P-0003698 chr13  32913797  32913797 G     C      1027   162 0.158 BRCA2  
#>  8 P-0003698 chr13  32914259  32914259 G     C      1021   182 0.178 BRCA2  
#>  9 P-0003698 chr19  11136104  11136104 G     T       573    98 0.171 SMARCA4
#> 10 P-0003698 chr22  41543840  41543840 G     A       416    45 0.108 EP300  
#> # ℹ 224,929 more rows
```
