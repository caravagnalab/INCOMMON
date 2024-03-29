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
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit

cox_fit = function(x, gene, tumor_type, survival_time, survival_status, covariates = c('age', 'sex', 'tmb')){
  
  data = prepare_km_fit_input(x, tumor_type, gene)
  data = data %>% dplyr::mutate(group = factor(class))
  data = data %>% dplyr::mutate(group = relevel(group, ref = grep('WT', unique(data$group), value = T)))
  
  
  formula = paste0('survival::Surv(',survival_time,', ',survival_status,') ~ group')
  
  for(c in covariates) {
    what = grep(c, colnames(data), ignore.case = T, value = TRUE)
    for(w in what) {
      if(is.numeric(data[[w]])){
        q =  quantile(data[w], na.rm = T)['50%']
        data[[w]] = ifelse(data[[w]] > q, paste0('>', round(q, 0)), paste0('<=', round(q, 0)))
        data[[w]] = factor(data[[w]])
        data[[w]] = relevel(data[[w]], ref = grep('<=', unique(data[[w]]), value = T))
      }
      formula = paste(formula, w, sep = ' + ')
    }
  }
  
  fit = survival::coxph(
    formula = formula %>% as.formula(),
    data = data %>%
      as.data.frame()
  )
  
  if(!('survival' %in% names(x))) x$survival = list()
  if(!(tumor_type %in% names(x$survival))) x$survival[[tumor_type]] = list()
  if(!(gene %in% names(x$survival[[tumor_type]]))) x$survival[[tumor_type]][[gene]] = list()
  
  x$survival[[tumor_type]][[gene]]$`cox_regression` = fit
  
  print(fit)
  
  return(x)
}
