#' Visualize frequency distribution of INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results, as
#' produced by using function `classify`.
#' @param tumor_type Tumor type for tumor-specific prior ('PANCA' for pan-cancer).
#' @param gene Gene for gene-specific prior.
#' @param ... Default S3 method parameter.
#' @return An object or a list of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Plot class fraction for a specific gene and tumour type
#' plot_class_fraction(x = MSK_classified, tumor_type = 'LUAD', gene = 'KRAS')
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom patchwork wrap_plots
plot_class_fraction = function(x, tumor_type, gene, ...){

  toplot = class_frequency(x, tumor_type, gene, ...)

  if('map_class' %in% colnames(toplot)){
    p = toplot %>%
      ggplot2::ggplot()+
      ggplot2::geom_bar(ggplot2::aes(x = '', y = frequency, fill = map_class), stat = 'identity')+
      scale_color_INCOMMON_class(aes = 'fill')+
      ggplot2::xlab('')+ggplot2::ylab('Fraction')+
      ggplot2::coord_flip()+
      INCOMMON:::my_ggplot_theme(cex = .8)+
      ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON state'))+
      ggplot2::labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))
  }

  if('relevant_class' %in% colnames(toplot)){
    p2 = toplot %>%
      dplyr::mutate(relevant_class = dplyr::case_when(
        grepl('with LOH', relevant_class) ~ 'with LOH',
        grepl('without LOH', relevant_class) ~ 'without LOH',
        grepl('with AMP', relevant_class) ~ 'with AMP',
        grepl('without AMP', relevant_class) ~ 'without AMP',
        grepl('Tier-2', relevant_class) ~ 'Tier-2',
      )) %>%
      ggplot2::ggplot()+
      ggplot2::geom_bar(ggplot2::aes(x = '', y = frequency, fill = relevant_class), stat = 'identity')+
      scale_color_INCOMMON_high_level_class(aes = 'fill')+
      ggplot2::xlab('')+ggplot2::ylab('Fraction')+
      ggplot2::coord_flip()+
      INCOMMON:::my_ggplot_theme(cex = .8)+
      ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON class'))+
      ggplot2::labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))

    p = patchwork::wrap_plots(p, p2, ncol = 1)
  }

  cli::cli_alert_info('The frequency of states {.field {toplot$state}} are {.field {round(toplot$frequency,2)}}')

  return(p)

}
