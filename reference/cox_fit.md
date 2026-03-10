# Fit multivariate Cox regression model based on INCOMMON classes.

Fit multivariate Cox regression model based on INCOMMON classes.

## Usage

``` r
cox_fit(
  x,
  gene,
  tumor_type,
  survival_time,
  survival_status,
  covariates = c("age", "sex", "tmb"),
  tmb_method = "median"
)
```

## Arguments

- x:

  An object of class `'INCOMMON'` containing the classification results
  as produced by function `classify`.

- gene:

  The gene on which patient's stratification is based.

- tumor_type:

  The tumor type of patients to stratify.

- survival_time:

  The variable in `clincal_data` to be used as survival time.

- survival_status:

  The variable in `clincal_data` to be used as survival status.

- covariates:

  The other covariates to be used in the mutlivariate regression.

- tmb_method:

  The method to define the reference value for tumor mutational burden
  TMB

## Value

An object of class `'INCOMMON'` containing an additional object
`survival`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Perform Cox regression based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
MSK_PAAD_output = cox_fit(x = MSK_PAAD_output, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS', covariates = c('age', 'sex', 'tmb'), tmb_method = ">10")
#> Joining with `by = join_by(id)`
#> [1] "Cox fit with INCOMMON groups:"
#> Call:
#> survival::coxph(formula = formula %>% stats::as.formula(), data = data %>% 
#>     as.data.frame())
#> 
#>                        coef exp(coef) se(coef)     z        p
#> groupBalanced Dosage 0.8179    2.2658   0.1159 7.058 1.69e-12
#> groupHigh Dosage     1.1569    3.1801   0.1242 9.318  < 2e-16
#> groupLow Dosage      0.6840    1.9818   0.1177 5.810 6.24e-09
#> 
#> Likelihood ratio test=102.1  on 3 df, p=< 2.2e-16
#> n= 1772, number of events= 1090 
#>    (7 observations deleted due to missingness)
#> [1] "Pairwise tests:"
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Fit: survival::coxph(formula = formula %>% stats::as.formula(), data = data %>% 
#>     as.data.frame())
#> 
#> Linear Hypotheses:
#>                                                  Estimate Std. Error z value
#> `groupHigh Dosage` - `groupBalanced Dosage` == 0  0.33900    0.08156   4.156
#> `groupLow Dosage` - `groupBalanced Dosage` == 0  -0.13391    0.07299  -1.835
#>                                                  Pr(>|z|)    
#> `groupHigh Dosage` - `groupBalanced Dosage` == 0 6.45e-05 ***
#> `groupLow Dosage` - `groupBalanced Dosage` == 0     0.123    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
```
