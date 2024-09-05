#' Visualize metastatic propnesity odds ratio in a volcano plot fashion.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results, as
#' produced by function `classify`.
#' @param tumor_type The tumour type for which classified data is available. If 'PANCA'
#' it pools from multiple tumour types into a pan-cancer visualisation.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Estimate the metastatic propensity associated with mutant TP53, PIK3CA and CDH1 with vs without CNA in BRCA.
#' for(g in c('TP53', 'PIK3CA', 'CDH1')){MSK_classified = met_propensity(x = MSK_classified, tumor_type = 'BRCA', gene = g)}
#' # Plot results in a volcano plot
#' plot_met_volcano(x = MSK_classified, tumor_type = 'BRCA')
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom ggrepel geom_label_repel

plot_met_volcano = function(x, tumor_type){
  stopifnot(inherits(x, 'INCOMMON'))
  stopifnot('metastatic_propensity' %in% names(x))

  if(tumor_type != 'PANCA'){
    stopifnot(tumor_type %in% names(x$metastatic_propensity))
    toplot = x$metastatic_propensity[[tumor_type]] %>%
      purrr::discard(is.null) %>%
      do.call(rbind, .)
  } else{
    toplot = lapply(x$metastatic_propensity, function(x) x) %>%
      unlist(recursive = F) %>%
      purrr::discard(is.null) %>%
      do.call(rbind, .)
  }

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
  ggplot2::ggplot(ggplot2::aes(
    x = log2(OR),
    y = -log10(p.value))) +
  ggplot2::geom_vline(xintercept = 1, linetype = 'longdash', color = 'indianred')+
  ggplot2::geom_vline(xintercept = -1, linetype = 'longdash', color = 'indianred')+
  ggplot2::geom_hline(yintercept = 1.3, linetype = 'longdash', color = 'indianred')+
  ggplot2::geom_point(ggplot2::aes(color = prevalent), shape = 19, size = 2)+
  ggrepel::geom_label_repel(data = toplot %>% dplyr::filter(prevalent != 'ns'),
                            ggplot2::aes(label = gene, fill = tumor_type),
                            size = 3)+
  ggplot2::xlim(-4,4)+
  scale_color_INCOMMON_high_level_class(aes = 'color')+
  scale_color_ttypes(aes = 'fill')+
  INCOMMON:::my_ggplot_theme(cex = .8)+
  ggplot2::xlab('Odds Ratio (log2)')+
  ggplot2::ylab('P-value (-log10)')+
  ggplot2::guides(fill = ggplot2::guide_legend(title = 'Tumor Type',
                                               override.aes = list(color=NA)
                                               ),
                  color = ggplot2::guide_legend(title = 'INCOMMON class'))
}
