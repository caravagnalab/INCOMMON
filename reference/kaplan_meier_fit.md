# Fit Kaplan-Meier survival model stratifying patients based on INCOMMON classes, for a single gene and tumour type.

Fit Kaplan-Meier survival model stratifying patients based on INCOMMON
classes, for a single gene and tumour type.

## Usage

``` r
kaplan_meier_fit(x, tumor_type, gene, survival_time, survival_status)
```

## Arguments

- x:

  An object of class `'INCOMMON'` containing the classification results
  as produced by function `classify`.

- tumor_type:

  The tumor type of patients to stratify.

- gene:

  The gene on which patient's stratification is based.

- survival_time:

  The variable in `clincal_data` to be used as survival time.

- survival_status:

  The variable in `clincal_data` to be used as survival status.

## Value

An object of class `'INCOMMON'` containing an additional object
`survival`.

## Examples

``` r
# First load example classified data
data(MSK_PAAD_output)
# Perform survival analysis based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
MSK_PAAD_output = kaplan_meier_fit(x = MSK_PAAD_output, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS')
#> Joining with `by = join_by(id)`
#> Call: survfit(formula = "survival::Surv(OS_MONTHS, OS_STATUS) ~ group", 
#>     data = data)
#> 
#>    7 observations deleted due to missingness 
#>                   n events median 0.95LCL 0.95UCL
#> WT              214     92   38.3   28.98    51.1
#> Low Dosage      541    344   17.4   15.34    19.8
#> Balanced Dosage 670    417   14.6   13.34    15.8
#> High Dosage     347    237   10.7    9.36    12.5
```
