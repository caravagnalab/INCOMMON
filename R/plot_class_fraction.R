#' Visualize frequency distribution of INCOMMON classes.
#'
#' @param x A list of objects of class \code{'INCOMMON'} containing the classification results, as 
#' produced by using function `classify`.
#' @return An object or a list of class \code{'ggplot2'}.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' data = data_MSK %>% filter(gene == 'KRAS', tumor_type 'PAAD')
#' samples = unique(data$sample)[1:10]
#' classified_data = lapply(samples, function(s){
#' x = init(mutations = data %>% filter(sample==s),
#'          sample = s,
#'          purity = data %>% filter(sample==s) %>% purity(),
#'          tumor_type = data %>% filter(sample==s) %>% tumor_type())
#' x = classify(
#'     x = x, 
#'     priors = NULL,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' })
#
#' plot_class_fraction(x, tumor_type = 'PAAD', gene = 'KRAS')
#' 
plot_class_fraction = function(x, tumor_type, gene){
  toplot = class_frequency(classified_data, tumor_type, gene)
  toplot %>% 
    ggplot()+
    geom_bar(aes(x = '', y = frequency, fill = state), stat = 'identity')+
    scale_color_INCOMMON_class(aes = 'fill')+
    xlab('')+ylab('Fraction')+
    coord_flip()+
    CNAqc:::my_ggplot_theme(cex = .8)+
    guides(fill = guide_legend(title = 'INCOMMON class'))+
    labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))
}