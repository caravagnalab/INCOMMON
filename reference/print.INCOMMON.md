# Print for class `'INCOMMON'`.

Print for class `'INCOMMON'`.

## Usage

``` r
# S3 method for class 'INCOMMON'
print(x, ...)
```

## Arguments

- x:

  An obj of class `'INCOMMON'`.

- ...:

  Default S3 method parameter.

## Value

Nothing.

## Examples

``` r
data(MSK_genomic_data)
data(MSK_clinical_data)
sample = 'P-0002081'
x = init(
  genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,],
  clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,]
)
#> ── INCOMMON - Inference of copy number and mutation multiplicity in oncology ───
#> 
#> ── Genomic data ──
#> 
#> ✔ Found 1 samples, with 4 mutations in 4 genes
#> 
#> ── Clinical data ──
#> 
#> ℹ Provided clinical features:
#> 
#> ✔ sample (required for classification)
#> ✔ purity (required for classification)
#> ✔ tumor_type
#> ✔ OS_MONTHS
#> ✔ OS_STATUS
#> ✔ SAMPLE_TYPE
#> ✔ MET_COUNT
#> ✔ METASTATIC_SITE
#> ✔ MET_SITE_COUNT
#> ✔ PRIMARY_SITE
#> ✔ SUBTYPE_ABBREVIATION
#> ✔ GENE_PANEL
#> ✔ SEX
#> ✔ TMB_NONSYNONYMOUS
#> ✔ FGA
#> ✔ AGE_AT_SEQUENCING
#> ✔ RACE
#> 
#> ✔ Found 1 matching samples
#> ✔ No mismatched samples
print(x)
#> ── [ INCOMMON ]  4 PASS mutations across 1 samples,
#> with 4 mutant genes across 1
#> ℹ Average sample purity: 0.6
#> ℹ Average sequencing depth: 380
#> # A tibble: 4 × 27
#>   sample    tumor_type purity chr      from     to ref   alt      DP    NV   VAF
#>   <chr>     <chr>       <dbl> <chr>   <dbl>  <dbl> <chr> <chr> <int> <int> <dbl>
#> 1 P-0002081 LUAD          0.6 chr12  2.54e7 2.54e7 C     A       743   378 0.509
#> 2 P-0002081 LUAD          0.6 chr17  7.58e6 7.58e6 G     A       246   116 0.472
#> 3 P-0002081 LUAD          0.6 chr19  1.22e6 1.22e6 C     A       260   122 0.469
#> 4 P-0002081 LUAD          0.6 chr19  1.11e7 1.11e7 -     C       271   133 0.491
#> # ℹ 16 more variables: gene <chr>, gene_role <chr>, OS_MONTHS <dbl>,
#> #   OS_STATUS <dbl>, SAMPLE_TYPE <chr>, MET_COUNT <dbl>, METASTATIC_SITE <chr>,
#> #   MET_SITE_COUNT <dbl>, PRIMARY_SITE <chr>, SUBTYPE_ABBREVIATION <chr>,
#> #   GENE_PANEL <chr>, SEX <chr>, TMB_NONSYNONYMOUS <dbl>, FGA <dbl>,
#> #   AGE_AT_SEQUENCING <dbl>, RACE <chr>
```
