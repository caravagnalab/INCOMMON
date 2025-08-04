#' Visualize metastatic propnesity odds ratio in a volcano plot fashion.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results, as
#' produced by function `classify`.
#' @param tumor_type The tumour type for which classified data is available.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Estimate the metastatic tropism associated with mutant TP53, PIK3CA and CDH1 with vs without CNA,
#' # in BRCA, with respect to Liver and Lymph
#' for(g in c('TP53', 'PIK3CA', 'CDH1')){for(m in c('Liver', 'Lymph')){MSK_PAAD_output = met_tropism(x = MSK_PAAD_output, gene = g, tumor_type = 'BRCA', metastatic_site = m)}}
#' # Plot results in a volcano plot
#' plot_tropism(x = MSK_PAAD_output, tumor_type = 'BRCA')
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom ggrepel geom_label_repel

plot_tropism = function(x, tumor_type){
  stopifnot(inherits(x, 'INCOMMON'))
  stopifnot('metastatic_tropism' %in% names(x))

  stopifnot(tumor_type %in% names(x$metastatic_tropism))
  toplot = x$metastatic_tropism[[tumor_type]] %>%
    unlist(recursive = F) %>%
    do.call(rbind, .)

  toplot = toplot %>%
    dplyr::mutate(
      prevalent = dplyr::case_when(
        OR >= 1 & p.value <= .05 & grepl('LOH', class) ~ 'with LOH',
        OR >= 1 & p.value <= .05 & grepl('AMP', class) ~ 'with AMP',
        OR < 1 & p.value <= .05 & grepl('LOH', class) ~ 'without LOH',
        OR < 1 & p.value <= .05 & grepl('AMP', class) ~ 'without LOH',
        TRUE ~ 'ns'
      )
    )

  toplot %>%
    dplyr::filter(!is.na(low),
                  !is.na(up)) %>%
    ggplot2::ggplot(ggplot2::aes(
      y = gene, x = log2(OR))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 'longdash', color = 'indianred')+
    ggplot2::geom_point(ggplot2::aes(color = prevalent), shape = 19, size = 2)+
    ggplot2::geom_errorbar(ggplot2::aes(xmin = log2(low),
                                        xmax = log2(up),
                                        color = prevalent), width = .25)+
    ggplot2::facet_wrap(~metastatic_site, scales = 'free_y')+
    scale_color_INCOMMON_high_level_class(aes = 'color')+
    my_ggplot_theme(cex = .8)+
    ggplot2::xlab('Odds Ratio (log2)')+
    ggplot2::ylab('')+
    ggplot2::guides(color = ggplot2::guide_legend(title = 'INCOMMON class'))
}
