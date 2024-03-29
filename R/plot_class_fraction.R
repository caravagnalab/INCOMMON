#' Visualize frequency distribution of INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results, as 
#' produced by using function `classify`.
#' @param tumor_type Tumor type for tumor-specific prior ('PANCA' for pan-cancer).
#' @param gene Gene for gene-specific prior.
#' @return An object or a list of class \code{'ggplot2'}.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
plot_class_fraction = function(x, tumor_type, gene){
  toplot = class_frequency(x, tumor_type, gene)
  toplot %>% 
    ggplot2::ggplot()+
    ggplot2::geom_bar(ggplot2::aes(x = '', y = frequency, fill = state), stat = 'identity')+
    scale_color_INCOMMON_class(aes = 'fill')+
    ggplot2::xlab('')+ggplot2::ylab('Fraction')+
    ggplot2::coord_flip()+
    CNAqc:::my_ggplot_theme(cex = .8)+
    ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON class'))+
    ggplot2::labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))
}