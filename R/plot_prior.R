#' Visualize prior distribution for a gene (tumor-specific or pancancer).
#'
#' @param x A prior distribution in the format required for \code{INCOMMON},
#' such as \code{INCOMMON::priors_k_m}.
#' @param tumor_type Tumor type for tumor-specific prior ('PANCA' for pan-cancer).
#' @param gene Gene for gene-specific prior.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Plot classification results for a specific sample
#' plot_prior(x = MSK_PAAD_output, gene = 'TP53', tumor_type = 'PAAD')
#' @importFrom dplyr filter mutate rename select %>%
plot_prior = function(x, gene, tumor_type){
  toplot = x %>%
    dplyr::filter(gene == !!gene, tumor_type == !!tumor_type)
  toplot %>%
    ggplot2::ggplot(ggplot2::aes(x = k, y = m,))+
    ggplot2::geom_tile(ggplot2::aes(fill = log10(n)))+
    ggplot2::scale_fill_viridis_c(option = 'turbo')+
    ggplot2::geom_text(ggplot2::aes(label = round(n, 0)), color = 'white')+
    ggplot2::labs(
      title = paste0(gene, ' (', tumor_type, ')'),
      subtitle = paste0('Total samples N = ', unique(round(toplot$N,0))))+
    my_ggplot_theme()
}
