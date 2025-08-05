# Multivariate Cox regression
#' Fit multivariate Cox regression model based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @param survival_time The variable in `clincal_data` to be used as survival time.
#' @param survival_status The variable in `clincal_data` to be used as survival status.
#' @param covariates The other covariates to be used in the mutlivariate regression.
#' @param tmb_method The method to define the reference value for tumor mutational burden TMB
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Perform Cox regression based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
#' MSK_PAAD_output = cox_fit(x = MSK_PAAD_output, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS', covariates = c('age', 'sex', 'tmb'), tmb_method = ">10")
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit
#' @importFrom stats relevel quantile as.formula

cox_fit = function(x, gene, tumor_type, survival_time, survival_status,
                   covariates = c('age', 'sex', 'tmb'),
                   tmb_method = 'median'){

  if(!("genotype" %in% (input(x) %>% names()))) x = mutant_dosage_classification(x)

  data = prepare_km_fit_input(x, tumor_type, gene)

  fit = function(data, baseline){
    if(baseline){
      # Baseline fit
      data = data %>%
        dplyr::mutate(class = dplyr::case_when(
          grepl('WT', class) ~ 'WT',
          TRUE ~ 'Mutant'
        ))
    } else {
      data = data %>%
        dplyr::mutate(class = dplyr::case_when(
          grepl('WT', class) ~ 'WT',
          TRUE ~ class
        ))
    }
    data = data %>% dplyr::mutate(group = factor(class))
    data = data %>% dplyr::mutate(group = relevel(group, ref = grep('WT', unique(data$group), value = T)))

    # Essential regression formula
    formula = paste0('survival::Surv(',survival_time,', ',survival_status,') ~ group')

    # Add covariates to regression formula
    for(c in covariates) {
      what = grep(c, colnames(data), fixed = T, value = TRUE)
      for(w in what) {
        if(grepl('tmb', w, ignore.case = T) & tmb_method != 'median'){
          data[[w]] = ifelse(data[[w]] > 10, '> 10', '<= 10')
          data[[w]] = factor(data[[w]])
          data[[w]] = stats::relevel(data[[w]], ref = '<= 10', value = T)
        }
        # if(is.numeric(data[[w]])){
        #   q =  stats::quantile(data[w], na.rm = T)['50%']
        #   data[[w]] = ifelse(data[[w]] > q, paste0('>', q), paste0('<=', q))
        #   data[[w]] = factor(data[[w]])
        #   data[[w]] = stats::relevel(data[[w]], ref = grep('<=', unique(data[[w]]), value = T))
        # }
        formula = paste(formula, w, sep = ' + ')
      }
    }

    fit = survival::coxph(
      formula = formula %>% stats::as.formula(),
      data = data %>%
        as.data.frame()
    )

    adj_cox_fit = fit
    adj_cox_fit$coefficients[is.na(adj_cox_fit$coefficients)] = 0

    pw_test_formula = case_when(
      baseline == TRUE ~ "groupMutant = 0",
      baseline == FALSE ~ c("`groupHigh Dosage` - `groupBalanced Dosage` = 0", "`groupLow Dosage` - `groupBalanced Dosage` = 0"),
    )

    pw_test = summary(multcomp::glht(adj_cox_fit, linfct = c(pw_test_formula)))

    ph_test = tryCatch(survival::cox.zph(adj_cox_fit), error = function(e) {return(tibble(NULL))})

    return(dplyr::tibble(cox_fit = list(adj_cox_fit), pw_test = list(pw_test), ph_test = list(ph_test), formula))
  }

  baseline_fit = fit(data = data, baseline = TRUE)
  incommon_fit = fit(data = data, baseline = FALSE)

  if(!('survival' %in% names(x))) x$survival = list()
  if(!(tumor_type %in% names(x$survival))) x$survival[[tumor_type]] = list()
  if(!(gene %in% names(x$survival[[tumor_type]]))) x$survival[[tumor_type]][[gene]] = list()

  x$survival[[tumor_type]][[gene]]$`cox_regression`$baseline = baseline_fit
  x$survival[[tumor_type]][[gene]]$`cox_regression`$incommon = incommon_fit

  print("Cox fit with INCOMMON groups:")
  print(incommon_fit$cox_fit[[1]])
  print("Pairwise tests:")
  print(incommon_fit$pw_test[[1]])

  return(x)
}
