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
#' data(MSK_PAAD_output)
#' # Plot class fraction for a specific gene and tumour type
#' MSK_PAAD_output = mutant_dosage_classification(MSK_PAAD_output)
#' plot_class_fraction(x = MSK_PAAD_output, tumor_type = 'PAAD', gene = 'KRAS')
#' @importFrom dplyr filter mutate rename select group_by reframe ungroup arrange desc %>%
#' @importFrom patchwork wrap_plots
plot_class_fraction = function(x, tumor_type = NULL, gene = NULL, ...){

  toplot = class_frequency(x, tumor_type, gene, ...)

  title = gene
  title = ifelse(!is.null(tumor_type), paste0(title, ' (', tumor_type, ')'), title)

    p = toplot %>%
      dplyr::mutate(
        class = factor(class, levels = c('Low Dosage', 'Balanced Dosage', 'High Dosage') %>% rev()),
        SAMPLE_TYPE = factor(SAMPLE_TYPE, levels = c('Primary', 'Metastasis'))) %>%
      dplyr::group_by(SAMPLE_TYPE) %>%
      dplyr::arrange(dplyr::desc(class), .by_group = T) %>%
      dplyr::mutate(z = cumsum(n)-.5*n) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot()+
      ggplot2::geom_bar(ggplot2::aes(x = '', y = n, fill = class), stat = 'identity')+
      ggplot2::geom_text(ggplot2::aes(x = '', y = z, label = paste0(100*round(f,2),'%')), color = 'white')+
      scale_color_INCOMMON_class(aes = 'fill')+
      ggplot2::xlab('')+ggplot2::ylab('N')+
      ggplot2::coord_flip()+
      my_ggplot_theme(cex = .8)+
      ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON state'))+
      ggplot2::labs(title = title, subtitle = paste0('Total Number of Samples N = ', unique(toplot$tot)))+
      ggplot2::facet_wrap(~SAMPLE_TYPE, scales = 'free')

  return(p)

}
