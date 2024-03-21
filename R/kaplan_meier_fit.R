#' Fit Kaplan-Meier survival model based on INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results for
#' multiple samples, as produced by using function `classify`.
#' @param tumor_type The selected tumor type.
#' @param gene The selected gene.
#' @return An table containing the selected gene, tumor_type, the data table used to
#' fit the model and an object of class \code{'survfit'}.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv survfit

kaplan_meier_fit = function(x, tumor_type, gene) {
  
  data = prepare_km_fit_input(x, tumor_type, gene)
  data = data %>%
    dplyr::mutate(group = factor(group,
                                 levels = c(
                                   grep('WT', unique(data$group), value = T),
                                   grep('Mutant', unique(data$group), value = T) %>% grep('with', ., invert = T, value = T),
                                   grep('Mutant', unique(data$group), value = T) %>% grep('with', ., value = T)
                                 )))
  fit = survival::survfit(formula = survival::Surv(OS_MONTHS, OS_STATUS) ~ group, data = data)
  
  return(tibble(gene = gene, tumor_type = tumor_type, data = list(data), fit = list(fit)))
}