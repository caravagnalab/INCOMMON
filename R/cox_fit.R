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
#' data(MSK_classified)
#' # Perform Cox regression based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
#' MSK_classified = cox_fit(x = MSK_classified, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS', covariates = c('age', 'sex', 'tmb'), tmb_method = ">10")
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit
#' @importFrom stats relevel quantile as.formula

cox_fit = function(x, gene, tumor_type, survival_time, survival_status,
                   covariates = c('age', 'sex', 'tmb'),
                   tmb_method = 'median'){

  if(!("genotype" %in% (classification(x) %>% names()))) x = genome_interpreter(x)

  data = prepare_km_fit_input(x, tumor_type, gene)

  fit = function(data, baseline){
    if(baseline){
      # Baseline fit
      data = data %>%
        dplyr::mutate(class = dplyr::case_when(
          grepl('WT', class) ~ class,
          TRUE ~ strsplit(class, split = ' with')[[1]][1]
        ))
    }
    data = data %>% dplyr::mutate(group = factor(class))
    data = data %>% dplyr::mutate(group = relevel(group, ref = grep('WT', unique(data$group), value = T)))

    # Essential regression formula
    formula = paste0('survival::Surv(',survival_time,', ',survival_status,') ~ group')

    # Add covariates to regression formula
    for(c in covariates) {
      what = grep(c, colnames(data), ignore.case = T, value = TRUE)
      for(w in what) {
        if(grepl('tmb', w, ignore.case = T) & tmb_method != 'median'){
          data[[w]] = ifelse(data[[w]] > 10, '> 10', '<= 10')
          data[[w]] = factor(data[[w]])
          data[[w]] = stats::relevel(data[[w]], ref = '<= 10', value = T)
        }
        if(is.numeric(data[[w]])){
          q =  stats::quantile(data[w], na.rm = T)['50%']
          data[[w]] = ifelse(data[[w]] > q, paste0('>', round(q, 0)), paste0('<=', round(q, 0)))
          data[[w]] = factor(data[[w]])
          data[[w]] = stats::relevel(data[[w]], ref = grep('<=', unique(data[[w]]), value = T))
        }
        formula = paste(formula, w, sep = ' + ')
      }
    }

    fit = survival::coxph(
      formula = formula %>% stats::as.formula(),
      data = data %>%
        as.data.frame()
    )

    return(fit)
  }

  baseline_fit = fit(data = data, baseline = TRUE)
  incommon_fit = fit(data = data, baseline = FALSE)

  if(!('survival' %in% names(x))) x$survival = list()
  if(!(tumor_type %in% names(x$survival))) x$survival[[tumor_type]] = list()
  if(!(gene %in% names(x$survival[[tumor_type]]))) x$survival[[tumor_type]][[gene]] = list()

  x$survival[[tumor_type]][[gene]]$`cox_regression`$baseline = baseline_fit
  x$survival[[tumor_type]][[gene]]$`cox_regression`$incommon = incommon_fit

  print(incommon_fit)

  return(x)
}
