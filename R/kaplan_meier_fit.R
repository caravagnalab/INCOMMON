#' Fit Kaplan-Meier survival model based on INCOMMON classes.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results
#' as produced by function `classify`.
#' @param tumor_type The tumor type of patients to stratify.
#' @param gene The gene on which patient's stratification is based.
#' @param survival_time The variable in `clincal_data` to be used as survival time.
#' @param survival_status The variable in `clincal_data` to be used as survival status.
#' @return An object of class \code{'INCOMMON'} containing an additional object `survival`.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit
#' @importFrom stats as.formula

kaplan_meier_fit = function(x, tumor_type, gene, survival_time, survival_status) {

  stopifnot(inherits(x, 'INCOMMON'))
  if(!("genotype" %in% (classification(x) %>% names()))) x = genome_interpreter(x)

  data = prepare_km_fit_input(x, tumor_type, gene)

  fit = function(data, baseline){
    if(baseline){
      data = data %>%
        dplyr::mutate(group = dplyr::case_when(
          grepl('WT', class) ~ class,
          TRUE ~ strsplit(class, split = ' with')[[1]][1]
        ))

      data = data %>%
        dplyr::mutate(group = factor(group,
                                     levels = c(
                                       grep('WT', unique(data$group), value = T),
                                       grep('Mutant', unique(data$group), value = T)
                                     )))
    } else{
      data = data %>%
        dplyr::mutate(group = factor(class,
                                     levels = c(
                                       grep('WT', unique(data$class), value = T),
                                       grep('Mutant', unique(data$class), value = T) %>% grep('without', ., value = T),
                                       grep('Mutant', unique(data$class), value = T) %>% grep('without', ., invert = T, value = T)
                                     )))
    }

    # Essential regression formula
    formula = paste0('survival::Surv(',survival_time,', ',survival_status,') ~ group')

    fit = survival::survfit(formula = formula %>% stats::as.formula(),
                            data = data)

    fit$call$formula = formula

    names(fit$strata) = gsub(fit$strata %>% names(),
                             pattern = 'group=',
                             replacement = '')

    return(fit)
  }

  baseline_fit = fit(data = data, baseline = TRUE)
  incommon_fit = fit(data = data, baseline = FALSE)

  if(!('survival' %in% names(x))) x$survival = list()
  if(!(tumor_type %in% names(x$survival))) x$survival[[tumor_type]] = list()
  if(!(gene %in% names(x$survival[[tumor_type]]))) x$survival[[tumor_type]][[gene]] = list()

  x$survival[[tumor_type]][[gene]]$`kaplan-meier` = list()

  x$survival[[tumor_type]][[gene]]$`kaplan-meier`$baseline = baseline_fit
  x$survival[[tumor_type]][[gene]]$`kaplan-meier`$incommon = incommon_fit

  print(incommon_fit)

  return(x)
}
