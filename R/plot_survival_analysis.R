#' Visualise survival analysis based on INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results for
#' multiple samples, as produced by using function `classify`.
#' @param tumor_type The selected tumour type.
#' @param gene The selected gene.
#' @param cox_covariates A character vector listing the covariates to be used in the multivariarte regression.
#' @return An object or a list of class \code{'ggplot2'} showing Kaplan-Meier curves and 
#' Cox regression forest plot.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom survminer ggsurvplot
#' @importFrom patchwork wrap_plots plot_layout

plot_survival_analysis = function(x, tumor_type, gene, cox_covariates = c('age', 'sex', 'tmb')){
  
  stopifnot(inherits(x, 'INCOMMON'))
  
  stopifnot('survival' %in% names(x))
  stopifnot(tumor_type %in% names(x$survival))
  stopifnot(gene %in% names(x$survival[[tumor_type]]))
  stopifnot('kaplan-meier' %in% names(x$survival[[tumor_type]][[gene]]))
  stopifnot('cox_regression' %in% names(x$survival[[tumor_type]][[gene]]))
  
  
  km_fit = x$survival[[tumor_type]][[gene]]$`kaplan-meier`
  
  km_data = prepare_km_fit_input(x, tumor_type, gene)
  
  km_data = km_data %>%
    dplyr::mutate(group = factor(class,
                                 levels = c(
                                   grep('WT', unique(km_data$class), value = T),
                                   grep('Mutant', unique(km_data$class), value = T) %>% grep('with', ., invert = T, value = T),
                                   grep('Mutant', unique(km_data$class), value = T) %>% grep('with', ., value = T)
                                 )))
  
  km_plot = survminer::ggsurvplot(
    fit = km_fit,
    censor = F,
    conf.int = F,
    data = km_data,
    ylab = "Overall Survival",
    xlab = "Time (months)",
    fontsize = 4,
    risk.table = TRUE,
    risk.table.col = "strata",
    table.fontsize = 0.1,
    ggtheme = CNAqc:::my_ggplot_theme(cex = .8),
    tables.theme = CNAqc:::my_ggplot_theme(cex = .8),
    palette = surv_colors(unique(km_data$gene_role))
  ) 
  
  km_plot$plot$data$tumor_type = unique(km_data$tumor_type)
  km_plot$data.survplot$tumor_type = unique(km_data$tumor_type)
  
  km_plot$plot = km_plot$plot + 
    ggplot2::xlab('') + 
    ggplot2::guides(color = 'none') + ggplot2::facet_wrap(~tumor_type)
  
  km_plot$table = km_plot$table + 
    ggplot2::ylab('') + 
    ggplot2::theme(legend.position = 'none')
  
cox_fit = x$survival[[tumor_type]][[gene]]$`cox_regression`

forest_plot = forest_plot(cox_fit)

patchwork::wrap_plots(km_plot$plot, 
                      km_plot$table, 
                      forest_plot)+
  patchwork::plot_layout(heights = c(2.5,1,2))
}
