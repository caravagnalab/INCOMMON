#' Visualise the posterior distribution on (k,m) configurations.
#' @param x An object of class INCOMMON.
#' @param k_max The maximum allowed value of total copyu number.
#' @param z_km The list of posterior distributions over (k,m) configurations.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Plot classification results for a specific sample
#' x = subset_sample(x = MSK_PAAD_output, sample_list = "P-0000142")
#' plot_posterior_k_m(x = x)
#' @importFrom dplyr group_by arrange desc slice_head
#' @importFrom ggh4x facet_nested_wrap
#' @importFrom ggpubr get_legend as_ggplot
plot_posterior_k_m = function(x, k_max = NULL, z_km = NULL){
  # x = subset_sample(x = x, sample_list = sample)
  inp = input(x)
  M = inp %>% nrow()

  if(is.null(k_max)){
    k_max = x$parameters$k_max
  }

  if(is.null(z_km)){
    z_km = inp$z_km
  }

  toplot = lapply(1:M, function(i){
    output = z_km[[i]]
    output$id = paste(inp[i,]$gene, inp[i,]$NV, sep = ":")
    output
  }) %>% do.call(rbind, .)

  toplot = toplot %>%
    tidyr::complete(
      k = 1:8,
      m = 1:8,
      id,
      fill = list(z_km = min(toplot$z_km))  # Or use a small dummy like 1e-40 if log-scaling
    )

  pg = toplot %>%
    ggplot2::ggplot(ggplot2::aes(x=k,y=m))+
    ggplot2::geom_tile(ggplot2::aes(fill = log10(z_km)))+
    ggplot2::scale_fill_viridis_c()+
    ggplot2::guides(fill = ggplot2::guide_colorbar(title = 'Posterior Probability (log10)', position = 'bottom'))+
    my_ggplot_theme()

  leg = ggpubr::get_legend(pg)
  leg = ggpubr::as_ggplot(leg)

  p = toplot %>%
    ggplot2::ggplot(ggplot2::aes(x=k,y=m))+
    ggplot2::geom_contour_filled(ggplot2::aes(z = log10(z_km)), bins = 50, alpha = 0.85)+
    ggplot2::guides(fill = 'none')+
    ggplot2::coord_cartesian(expand = FALSE, clip = 'off')+
    ggplot2::geom_point(
      data = toplot %>%
        dplyr::group_by(id) %>% arrange(desc(z_km), .by_group = T) %>% slice_head(n=1),
      color = 'firebrick',
      shape = 8)+
    ggplot2::scale_y_continuous(position = 'right') +
    ggplot2::facet_wrap(~id, ncol = 1, strip.position = 'left')+
    my_ggplot_theme()+
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_blank(),panel.grid = ggplot2::element_line(colour = 'black'),
      legend.position = 'bottom',
      axis.text = ggplot2::element_text(),
      axis.title = ggplot2::element_text()
    )+
    ggplot2::labs(x = 'Total CN (k)', y = 'Multiplicity (m)')

  p = p/leg+plot_layout(heights = c(10,1))

  p
}
