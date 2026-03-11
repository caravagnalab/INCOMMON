# Group patients by gene mutant mutant dosage using gene-role specific thresholds.

Group patients by gene mutant mutant dosage using gene-role specific
thresholds.

## Usage

``` r
mutant_dosage_classification(
  x,
  TSG_low = 0.25,
  TSG_high = 0.75,
  ONC_low = 0.33,
  ONC_high = 0.66
)
```

## Arguments

- x:

  An object of class INCOMMON.

- TSG_low:

  The lower cutoff for mutant dosage classification of tumour suppressor
  genes.

- TSG_high:

  The upper cutoff for mutant dosage classification of tumour suppressor
  genes.

- ONC_low:

  The lower cutoff for mutant dosage classification for oncogenes.

- ONC_high:

  The upper cutoff for mutant dosage classification for oncogenes.

## Value

An object of class INCOMMON with new columns reporting mean FAM and
class assignment.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
mutant_dosage_classification(MSK_PAAD_output, TSG_low = .25, TSG_high = .75, ONC_low = .33, ONC_high = .66)
#> Joining with `by = join_by(id)`
#> ── [ INCOMMON ]  7839 PASS mutations across 1779 samples,
#> with 280 mutant genes 
#> ℹ Average sample purity: 0.4
#> ℹ Average sequencing depth: 649
#> # A tibble: 7,839 × 192
#>    sample genotype tumor_type purity chr     from     to ref   alt      NV    DP
#>    <chr>  <chr>    <chr>       <dbl> <chr>  <dbl>  <dbl> <chr> <chr> <int> <int>
#>  1 P-000… KRAS wi… PAAD          0.4 chr12 2.54e7 2.54e7 C     C       273  1404
#>  2 P-000… KRAS wi… PAAD          0.4 chr17 7.58e6 7.58e6 G     G        53   671
#>  3 P-000… KRAS wi… PAAD          0.4 chr2  4.77e7 4.77e7 T     T        31   481
#>  4 P-000… KRAS wi… PAAD          0.4 chr5  1.28e6 1.28e6 G     G        34   227
#>  5 P-000… KRAS wi… PAAD          0.8 chr12 2.54e7 2.54e7 C     C       474   941
#>  6 P-000… KRAS wi… PAAD          0.8 chr5  1.12e8 1.12e8 G     G       164   424
#>  7 P-000… KRAS wi… PAAD          0.8 chr11 8.60e7 8.60e7 T     T       210   601
#>  8 P-000… KRAS wi… PAAD          0.8 chr13 3.29e7 3.29e7 TC    TC      160   493
#>  9 P-000… BRAF wi… PAAD          0.6 chr7  1.40e8 1.40e8 A     A       308   736
#> 10 P-000… BRAF wi… PAAD          0.6 chr1  1.15e8 1.15e8 T     T       188   506
#> # ℹ 7,829 more rows
#> # ℹ 181 more variables: gene <chr>, HGVSp_Short <chr>, Entrez_Gene_Id <dbl>,
#> #   Center <chr>, NCBI_Build <chr>, Chromosome <chr>, Strand <chr>,
#> #   Consequence <chr>, Variant_Classification <chr>, Variant_Type <chr>,
#> #   Tumor_Seq_Allele2 <chr>, dbSNP_RS <chr>, dbSNP_Val_Status <lgl>,
#> #   Matched_Norm_Sample_Barcode <lgl>, Match_Norm_Seq_Allele1 <lgl>,
#> #   Match_Norm_Seq_Allele2 <lgl>, Tumor_Validation_Allele1 <lgl>, …
```
