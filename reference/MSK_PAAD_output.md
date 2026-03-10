# Data from the MSK-MetTropism cohort classified with INCOMMON

This is an INCOMMON object resulting from the classification of all the
example data.

## Usage

``` r
data(MSK_PAAD_output)
```

## Format

An INCOMMON object with genomics, clinical data and classification
results.

## Source

"msk_classified_with_priors.rds" at Zenodo
https://zenodo.org/records/10927218

## Examples

``` r
data(MSK_PAAD_output)
MSK_PAAD_output
#> ── [ INCOMMON ]  175054 PASS mutations across 24018 samples, with 290 mutant gen
#> ℹ Average sample purity: 0.4
#> ℹ Average sequencing depth: 649
#> # A tibble: 7,839 × 186
#>    sample    tumor_type purity chr     from     to ref   alt      NV    DP gene 
#>    <chr>     <chr>       <dbl> <chr>  <dbl>  <dbl> <chr> <chr> <int> <int> <chr>
#>  1 P-0000142 PAAD          0.4 chr12 2.54e7 2.54e7 C     C       273  1404 KRAS 
#>  2 P-0000142 PAAD          0.4 chr17 7.58e6 7.58e6 G     G        53   671 TP53 
#>  3 P-0000142 PAAD          0.4 chr2  4.77e7 4.77e7 T     T        31   481 MSH2 
#>  4 P-0000142 PAAD          0.4 chr5  1.28e6 1.28e6 G     G        34   227 TERT 
#>  5 P-0000783 PAAD          0.8 chr12 2.54e7 2.54e7 C     C       474   941 KRAS 
#>  6 P-0000783 PAAD          0.8 chr5  1.12e8 1.12e8 G     G       164   424 APC  
#>  7 P-0000783 PAAD          0.8 chr11 8.60e7 8.60e7 T     T       210   601 EED  
#>  8 P-0000783 PAAD          0.8 chr13 3.29e7 3.29e7 TC    TC      160   493 BRCA2
#>  9 P-0000879 PAAD          0.6 chr7  1.40e8 1.40e8 A     A       308   736 BRAF 
#> 10 P-0000879 PAAD          0.6 chr1  1.15e8 1.15e8 T     T       188   506 NRAS 
#> # ℹ 7,829 more rows
#> # ℹ 175 more variables: HGVSp_Short <chr>, Entrez_Gene_Id <dbl>, Center <chr>,
#> #   NCBI_Build <chr>, Chromosome <chr>, Strand <chr>, Consequence <chr>,
#> #   Variant_Classification <chr>, Variant_Type <chr>, Tumor_Seq_Allele2 <chr>,
#> #   dbSNP_RS <chr>, dbSNP_Val_Status <lgl>, Matched_Norm_Sample_Barcode <lgl>,
#> #   Match_Norm_Seq_Allele1 <lgl>, Match_Norm_Seq_Allele2 <lgl>,
#> #   Tumor_Validation_Allele1 <lgl>, Tumor_Validation_Allele2 <lgl>, …
```
