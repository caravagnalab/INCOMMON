# Multivariate Cox regression

#' Fit multivariate Cox regression model based on INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results for
#' multiple samples, as produced by using function `classify`.
#' @param tumor_type The selected tumor type.
#' @param gene The selected gene.
#' @param covariates Covariates used in the multivariate regression.
#' @return An object of class \code{coxph}.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survival Surv coxph

cox_fit = function(x, gene, tumor_type, covariates = c('age', 'sex', 'tmb')){
  
  x = prepare_km_fit_input(x, tumor_type, gene)
  
  formula = 'survival::Surv(OS_MONTHS, OS_STATUS) ~ group'
  
  for(c in covariates) {
    what = grep(c, colnames(x), ignore.case = T, value = TRUE)
    for(w in what) {
      if(is.numeric(x[[w]])){
        q =  quantile(x[w], na.rm = T)['50%']
        x[[w]] = ifelse(x[[w]] > q, paste0('>', round(q, 0)), paste0('<=', round(q, 0)))
        x[[w]] = factor(x[[w]])
        x[[w]] = relevel(x[[w]], ref = grep('<=', unique(x[[w]]), value = T))
      }
      formula = paste(formula, w, sep = ' + ')
    }
  }
  
  fit = survival::coxph(
    formula = formula %>% as.formula(),
    data = x %>%
      dplyr::mutate(group = factor(group)) %>% 
      dplyr::mutate(group = relevel(group, ref = grep('WT', unique(x$group), value = T))) %>% 
      as.data.frame()
  )
  
  return(fit)
}
