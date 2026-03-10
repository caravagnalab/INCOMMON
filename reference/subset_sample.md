# Subset an INCOMMON object by sample ID.

Subset an INCOMMON object by sample ID.

## Usage

``` r
subset_sample(x, sample_list)
```

## Arguments

- x:

  An object of class `'INCOMMON'` generated with function `init`.

- sample_list:

  a list of identifiers for the samples to be subsetted

## Value

An object of class `INCOMMON` containing a subset of the original input.

## Examples

``` r
# First load example data
data(MSK_PAAD_output)
x = subset_sample(x = MSK_PAAD_output, sample_list = c("P-0000142"))
print(x)
#> ── [ INCOMMON ]  175054 PASS mutations across 24018 samples, with 290 mutant gen
#> ℹ Average sample purity: 0.4
#> ℹ Average sequencing depth: 649
#> # A tibble: 4 × 186
#>   sample    tumor_type purity chr      from     to ref   alt      NV    DP gene 
#>   <chr>     <chr>       <dbl> <chr>   <dbl>  <dbl> <chr> <chr> <int> <int> <chr>
#> 1 P-0000142 PAAD          0.4 chr12  2.54e7 2.54e7 C     C       273  1404 KRAS 
#> 2 P-0000142 PAAD          0.4 chr17  7.58e6 7.58e6 G     G        53   671 TP53 
#> 3 P-0000142 PAAD          0.4 chr2   4.77e7 4.77e7 T     T        31   481 MSH2 
#> 4 P-0000142 PAAD          0.4 chr5   1.28e6 1.28e6 G     G        34   227 TERT 
#> # ℹ 175 more variables: HGVSp_Short <chr>, Entrez_Gene_Id <dbl>, Center <chr>,
#> #   NCBI_Build <chr>, Chromosome <chr>, Strand <chr>, Consequence <chr>,
#> #   Variant_Classification <chr>, Variant_Type <chr>, Tumor_Seq_Allele2 <chr>,
#> #   dbSNP_RS <chr>, dbSNP_Val_Status <lgl>, Matched_Norm_Sample_Barcode <lgl>,
#> #   Match_Norm_Seq_Allele1 <lgl>, Match_Norm_Seq_Allele2 <lgl>,
#> #   Tumor_Validation_Allele1 <lgl>, Tumor_Validation_Allele2 <lgl>,
#> #   Match_Norm_Validation_Allele1 <lgl>, Match_Norm_Validation_Allele2 <lgl>, …
```
