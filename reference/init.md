# Prepare input for analyses with `'INCOMMON'`.

Process input data into an object of class `'INCOMMON'`, ready for
downstream analyses (e.g. `classify`).

## Usage

``` r
init(genomic_data, clinical_data, gene_roles = INCOMMON::cancer_gene_census)
```

## Arguments

- genomic_data:

  a data table of annotated mutations with columns sample name `sample`,
  mutant chromosome `chr`, mutation start position `from`, mutation end
  position `to`, reference allele `ref`, alternative allele `alt`,
  integer sequencing depth `DP`, integer number of reads with variant
  `NV`, variant allele frequency `VAF` and gene name `gene` as Hugo
  Symbol, protein sequence of the variant in HGVS recommended format,
  preferably 1-letter amino-acid code `HGVSp_Short`.

- clinical_data:

  a data table of clinical data with compulsory matching sample names
  `sample` and sample purity `purity`, and optional clinical features
  like tumor type ONCOTREE code `tumor_type` (required for tumor
  specific priors), overall survival status `OS_STATUS` and time
  `OS_MONTHS` (required for survival analysis), `SAMPLE_TYPE` (Primary
  or Metastasis) and number of metastases `MET_COUNT` (required for
  metastatic propensity analysis), metastatic site `METASTATIC_SITE`
  (required for metastatic tropism analysis), plus any other useful
  covariate.

- gene_roles:

  A data table reporting `gene` names and associated `gene_role`
  ("oncogene" or "TSG"). The default is taken from COSMIC Cancer Gene
  Census v98.

## Value

An object of class `INCOMMON`.

## Examples

``` r
# Example input data from the MSK-MetTropism cohort, released with the package
data(MSK_genomic_data)
print(MSK_genomic_data)
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
data(MSK_clinical_data)
print(MSK_clinical_data)
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
# Initialize the INCOMMON object (note the outputs to screen)
x = init(genomic_data = MSK_genomic_data, clinical_data = MSK_clinical_data)
#> ── INCOMMON - Inference of copy number and mutation multiplicity in oncology ───
#> 
#> ── Genomic data ──
#> 
#> ✔ Found 25659 samples, with 224939 mutations in 491 genes
#> ! No read counts found for 1393 mutations in 1393 samples
#> ! Gene name not provided for 1393 mutations
#> ! 201 genes could not be assigned a role (TSG or oncogene)
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
#> ✔ Found 25257 matching samples
#> ✖ Found 513 unmatched samples
# An S3 method can be used to report to screen what is in the object
print(x)
#> ── [ INCOMMON ]  223546 PASS mutations across 24266 samples, with 490 mutant gen
#> ℹ Average sample purity: 0.4
#> ℹ Average sequencing depth: 660
#> # A tibble: 223,546 × 27
#>    sample    tumor_type purity chr     from     to ref   alt      DP    NV   VAF
#>    <chr>     <chr>       <dbl> <chr>  <dbl>  <dbl> <chr> <chr> <int> <int> <dbl>
#>  1 P-0028912 CHOL          0.3 chr17 7.58e6 7.58e6 G     A       837   133 0.159
#>  2 P-0028912 CHOL          0.3 chr6  1.12e8 1.12e8 -     A       698   141 0.202
#>  3 P-0028912 CHOL          0.3 chrX  5.32e7 5.32e7 G     A       832    85 0.102
#>  4 P-0003698 BLCA          0.2 chr17 7.58e6 7.58e6 C     A       437   109 0.249
#>  5 P-0003698 BLCA          0.2 chr3  4.99e7 4.99e7 C     A       591    86 0.146
#>  6 P-0003698 BLCA          0.2 chr5  1.49e8 1.49e8 C     T       360    36 0.1  
#>  7 P-0003698 BLCA          0.2 chr13 3.29e7 3.29e7 G     C      1027   162 0.158
#>  8 P-0003698 BLCA          0.2 chr13 3.29e7 3.29e7 G     C      1021   182 0.178
#>  9 P-0003698 BLCA          0.2 chr19 1.11e7 1.11e7 G     T       573    98 0.171
#> 10 P-0003698 BLCA          0.2 chr22 4.15e7 4.15e7 G     A       416    45 0.108
#> # ℹ 223,536 more rows
#> # ℹ 16 more variables: gene <chr>, gene_role <chr>, OS_MONTHS <dbl>,
#> #   OS_STATUS <dbl>, SAMPLE_TYPE <chr>, MET_COUNT <dbl>, METASTATIC_SITE <chr>,
#> #   MET_SITE_COUNT <dbl>, PRIMARY_SITE <chr>, SUBTYPE_ABBREVIATION <chr>,
#> #   GENE_PANEL <chr>, SEX <chr>, TMB_NONSYNONYMOUS <dbl>, FGA <dbl>,
#> #   AGE_AT_SEQUENCING <dbl>, RACE <chr>
```
