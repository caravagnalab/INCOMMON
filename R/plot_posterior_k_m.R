#' Visualise the posterior distribution on (k,m) configurations.
#' @param x An object of class INCOMMON.
#' @param k_max The maximum allowed value of total copyu number.
#' @param z_km The list of posterior distributions over (k,m) configurations.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_genomic_data)
#' data(MSK_clinical_data)
#' data(cancer_gene_census)
#' data(priors_eta)
#' data(priors_k_m)
#' x = init(genomic_data = MSK_genomic_data, clinical_data = MSK_clinical_data,  gene_roles = cancer_gene_census)
#' x = subset_sample(x, "P-0028912")
#' x = classify(x = x, num_cores = 1, iter_warmup = 10, iter_sampling = 20, num_chains = 1)
#' # Plot classification results for a specific sample
#' plot_posterior_k_m(x = x, sample = "P-0028912", purity_error = 0.05)
#' @importFrom dplyr filter mutate rename select %>% tibble
plot_posterior_k_m = function(x, k_max = NULL, z_km = NULL){
  # x = subset_sample(x = x, sample_list = sample)
  inp = input(x)
  M = inp %>% nrow()

  if(is.null(k_max)){
    k_max = x$parameters$k_max
  }

  if(is.null(z_km)){
    z_km = x$input$z_km
  }

  toplot = lapply(1:M, function(i){
    output = z_km[[i]]
    output$id = paste(inp[i,]$gene, inp[i,]$NV, sep = ":")
    output
  }) %>% do.call(rbind, .)

  toplot %>% ggplot2::ggplot(ggplot2::aes(
    x = factor(k),
    y = factor(m),
    fill = log10(z_km)
  )) +
    ggplot2::geom_tile() +
    ggplot2::geom_point(
      data = toplot %>% dplyr::group_by(id) %>% dplyr::arrange(dplyr::desc(z_km)) %>% dplyr::slice_head(n=1),
      fill = 'firebrick',
      shape = 21,
      stroke = 0
    )+
    ggplot2::scale_fill_viridis_c(labels = scales::math_format(10 ^ .x)) +
    ggplot2::scale_x_discrete(breaks = scales::pretty_breaks(n=3))+
    ggplot2::scale_y_discrete(breaks = scales::pretty_breaks(n=3))+
    ggh4x::facet_nested_wrap( ~ id, ncol = 1, strip.position = 'left') +
    my_ggplot_theme() +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, margin = ggplot2::margin()),
    )+
    ggplot2::labs(
      x = "Total CN (k)", y = "Multiplicity (m)", fill = "Posterior Probability (log10)") +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = ggplot2::unit(2.5, 'cm')))
}
