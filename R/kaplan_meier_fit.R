#' Fit Kaplan-Meier survival model based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit

kaplan_meier_fit = function(x, tumor_type, gene) {
  
  stopifnot(inherits(x, 'INCOMMON'))
  if(!("genotype" %in% (classification(x) %>% names()))) x = genome_interpreter(x)
  
  data = prepare_km_fit_input(x, tumor_type, gene)
  
  data = data %>%
    dplyr::mutate(group = factor(class,
                                 levels = c(
                                   grep('WT', unique(data$class), value = T),
                                   grep('Mutant', unique(data$class), value = T) %>% grep('with', ., invert = T, value = T),
                                   grep('Mutant', unique(data$class), value = T) %>% grep('with', ., value = T)
                                 )))
  
  fit = survival::survfit(formula = survival::Surv(OS_MONTHS, OS_STATUS) ~ group, data = data)
  
  names(fit$strata) = gsub(fit$strata %>% names(), pattern='group=',replacement='')
  
  if(!('survival' %in% names(x))) x$survival = list()
  if(!(tumor_type %in% names(x$survival))) x$survival[[tumor_type]] = list()
  if(!(gene %in% names(x$survival[[tumor_type]]))) x$survival[[tumor_type]][[gene]] = list()
  
  x$survival[[tumor_type]][[gene]]$`kaplan-meier` = fit
  
  print(fit)
  
  return(x)
}