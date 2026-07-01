# Classify mutations using a Beta-Binomial model-based test.

Classify mutations using a Beta-Binomial model-based test.

## Usage

``` r
classify(
  x,
  k_max = 8,
  priors_k_m = priors_pcawg_hmf,
  priors_eta = priors_eta,
  purity_error = 0.05,
  num_cores = NULL,
  iter_warmup = 500,
  iter_sampling = 1000,
  num_chains = 4
)
```

## Arguments

- x:

  An object of class `'INCOMMON'` generated with function `init`.

- k_max:

  The maximum value of total copy number to be included in the model.

- priors_k_m:

  A dplyr::tibble or data frame with columns `gene`, `tumor_type`, `k`,
  `m`, `N` and `n` and `p` indicating tumor-specific or pan-cancer
  (PANCA) prior probabilities.

- priors_eta:

  A dplyr::tibble or data frame with columns `tumor_type`,`mean_eta`,
  `var_eta`, `N`,`alpha_eta` and `beta_eta` providing parameters of the
  Gamma prior distribution over the per copy sequencing rate.

- purity_error:

  The expected error on the input sample purity estimate.

- num_cores:

  The number of cores to use for parallel stan sampling.

- iter_warmup:

  The number of iterations of the stan warmup phase.

- iter_sampling:

  The number of iterations of the stan sampling phase.

- num_chains:

  The number of MCMC chains to be run in parallel.

## Value

An object of class `INCOMMON` containing the original input plus the
classification data and parameters.

## Examples

``` r
# First load example data
data(MSK_genomic_data)
data(MSK_clinical_data)
data(priors_pcawg_hmf)
data(priors_eta)
# Initialize the INCOMMON object for a single sample (note the outputs to screen)
sample = 'P-0002081'
x = init(genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,], clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,])
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
# Run INCOMMON classification
x = classify(x = x, priors_k_m = priors_pcawg_hmf, priors_eta = priors_eta, num_cores = 1, iter_warmup = 10, iter_sampling = 10, num_chains = 1)
#> ℹ Performing classification
#> 
#> ── INCOMMON inference of copy number and mutation multiplicity for sample P-0002
#> 
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 30, column 18 to column 33:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the division as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 34, column 19 to column 34:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the division as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 51, column 11 to column 26:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the division as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 65, column 18 to column 33:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the division as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 67, column 18 to column 33:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the division as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Warning in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 84, column 11 to column 26:
#>     Found int division:
#>         k_max * (k_max + 1) / 2
#>     Values will be rounded towards zero. If rounding is not desired you can
#>     write the di
#> vision as
#>         k_max * (k_max + 1) / 2.0
#>     If rounding is intended please use the integer division operator %/%.
#> Running MCMC with 1 chain...
#> 
#> Chain 1 WARNING: No variance estimation is 
#> Chain 1          performed for num_warmup < 20 
#> Chain 1 Iteration:  1 / 20 [  5%]  (Warmup) 
#> Chain 1 Iteration: 11 / 20 [ 55%]  (Sampling) 
#> Chain 1 Iteration: 20 / 20 [100%]  (Sampling) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is inf, but must be positive finite! (in '/tmp/Rtmp9dgjxH/model-1c86210cf57d.stan', line 43, column 2 to column 29)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 finished in 0.1 seconds.
# An S3 method can be used to report to screen what is in the object
print(x)
#> ── [ INCOMMON ]  4 PASS mutations across 1 samples,
#> with 4 mutant genes across 1
#> ℹ Average sample purity: 0.6
#> ℹ Average sequencing depth: 380
#> # A tibble: 4 × 36
#>   sample    tumor_type purity purity_map eta_map chr     from     to gene  ref  
#>   <chr>     <chr>       <dbl>      <dbl>   <dbl> <chr>  <dbl>  <dbl> <chr> <chr>
#> 1 P-0002081 LUAD          0.6      0.796    127. chr12 2.54e7 2.54e7 KRAS  C    
#> 2 P-0002081 LUAD          0.6      0.796    127. chr17 7.58e6 7.58e6 TP53  G    
#> 3 P-0002081 LUAD          0.6      0.796    127. chr19 1.22e6 1.22e6 STK11 C    
#> 4 P-0002081 LUAD          0.6      0.796    127. chr19 1.11e7 1.11e7 SMAR… -    
#> # ℹ 26 more variables: alt <chr>, NV <int>, DP <int>, VAF <dbl>, map_k <int>,
#> #   map_m <int>, gene_role <chr>, OS_MONTHS <dbl>, OS_STATUS <dbl>,
#> #   SAMPLE_TYPE <chr>, MET_COUNT <dbl>, METASTATIC_SITE <chr>,
#> #   MET_SITE_COUNT <dbl>, PRIMARY_SITE <chr>, SUBTYPE_ABBREVIATION <chr>,
#> #   GENE_PANEL <chr>, SEX <chr>, TMB_NONSYNONYMOUS <dbl>, FGA <dbl>,
#> #   AGE_AT_SEQUENCING <dbl>, RACE <chr>, bayes_p_purity <dbl>,
#> #   bayes_p_eta <dbl>, post_pred_p.value_DP <dbl>, …
```
