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
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Perform survival analysis based on the classification of KRAS mutant samples of pancreatic adenocarcinoma
#' MSK_classified = kaplan_meier_fit(x = MSK_classified, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS')
#' # Perform Cox regression
#' MSK_classified = cox_fit(x = MSK_classified, tumor_type = 'PAAD', gene = 'KRAS', survival_time = 'OS_MONTHS', survival_status = 'OS_STATUS', covariates = c('age', 'sex', 'tmb'), tmb_method = ">10")
#' plot_survival_analysis(x = MSK_classified, tumor_type = 'PAAD', gene = 'KRAS')
#' @importFrom dplyr filter mutate rename select %>% full_join
#' @importFrom survminer ggsurvplot
#' @importFrom patchwork wrap_plots plot_layout

plot_survival_analysis = function(x, tumor_type, gene, cox_covariates = c('age', 'sex', 'tmb')){

  stopifnot(inherits(x, 'INCOMMON'))

  stopifnot('survival' %in% names(x))
  stopifnot(tumor_type %in% names(x$survival))
  stopifnot(gene %in% names(x$survival[[tumor_type]]))
  stopifnot('kaplan-meier' %in% names(x$survival[[tumor_type]][[gene]]))
  stopifnot('cox_regression' %in% names(x$survival[[tumor_type]][[gene]]))

  km_data = prepare_km_fit_input(x, tumor_type, gene)
  baseline_km_fit = x$survival[[tumor_type]][[gene]]$`kaplan-meier`$baseline
  incommon_km_fit = x$survival[[tumor_type]][[gene]]$`kaplan-meier`$incommon

  plot = function(data, fit, baseline){
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
    } else {
      data = data %>%
        dplyr::mutate(group = factor(class,
                                     levels = c(
                                       grep('WT', unique(data$class), value = T),
                                       grep('Mutant', unique(data$class), value = T) %>% grep('without', ., value = T),
                                       grep('Mutant', unique(data$class), value = T) %>% grep('without', ., invert = T, value = T)
                                     )))
    }

    plot = survminer::ggsurvplot(
      fit = fit,
      censor = F,
      conf.int = F,
      data = data,
      ylab = "Overall Survival",
      xlab = "Time (months)",
      fontsize = 4,
      risk.table = TRUE,
      risk.table.col = "strata",
      risk.table.fontsize = 3,
      ggtheme = my_ggplot_theme(cex = .8),
      tables.theme = my_ggplot_theme(cex = .8),
      risk.table.y.text = FALSE,
      palette = surv_colors(unique(km_data$gene_role), baseline)
    )

    plot$plot$data$id = paste0(unique(km_data$gene), " (",unique(km_data$tumor_type),")")
    plot$data.survplot$id = paste0(unique(km_data$gene), " (",unique(km_data$tumor_type),")")

    plot$plot = plot$plot +
      ggplot2::xlab('') +
      ggplot2::guides(color = 'none') +
      ggplot2::facet_wrap(~id, strip.position = 'right')

    plot$table = plot$table +
      ggplot2::ylab('') +
      ggplot2::theme(legend.position = 'none', title = ggplot2::element_blank())

    return(plot)

  }

  # Baseline KM plot
  baseline_km_plot = plot(data = km_data,
                       fit = baseline_km_fit,
                       baseline = TRUE)
  # INCOMMON KM plot
  incommon_km_plot = plot(data = km_data,
                       fit = incommon_km_fit,
                       baseline = FALSE)

  baseline_cox_fit = x$survival[[tumor_type]][[gene]]$`cox_regression`$baseline
  incommon_cox_fit = x$survival[[tumor_type]][[gene]]$`cox_regression`$incommon

  baseline_forest_plot = forest_plot(baseline_cox_fit, baseline = TRUE)
  incommon_forest_plot = forest_plot(incommon_cox_fit, baseline = FALSE)

  baseline_forest_plot$data = baseline_forest_plot$data %>%
    full_join(incommon_forest_plot$data %>%
                select(var) %>%
                unique())

  incommon_levels = incommon_forest_plot$data$var %>% levels()

  incommon_forest_plot$data = incommon_forest_plot$data %>%
    full_join(baseline_forest_plot$data %>%
                select(var) %>%
                unique())

  levels = baseline_forest_plot$data$var %>% levels()


  levels[(length(levels)-2):length(levels)] =
    incommon_levels[(length(incommon_levels)-2):length(incommon_levels)]

  baseline_forest_plot$data$var = factor(baseline_forest_plot$data$var, levels = levels)
  incommon_forest_plot$data$var = factor(incommon_forest_plot$data$var, levels = levels)

  # Remove unnecessary labels
  incommon_forest_plot = incommon_forest_plot+
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())

  incommon_km_plot$plot = incommon_km_plot$plot+
    ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())

  patchwork::wrap_plots(baseline_km_plot$plot,
                        incommon_km_plot$plot,
                        baseline_km_plot$table,
                        incommon_km_plot$table,
                        baseline_forest_plot,
                        incommon_forest_plot,
                        design = 'AB\nAB\nAB\nCD\nEF\nEF')
}
