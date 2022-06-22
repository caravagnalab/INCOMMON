plot_fit = function(fit, target_gene) {
  # Get class of mutations on target gene
  gene_status = fit %>%
    dplyr::filter(gene == target_gene) %>%
    dplyr::mutate(class = paste0(class, ' (', nv, '/', dp, ')')) %>%
    dplyr::pull(class) %>%
    paste(collapse = ', ')
  
  # Get VAF of mutations on target gene
  gene_vafs = fit %>%
    dplyr::filter(gene == target_gene) %>%
    dplyr::pull(VAF)
  
  # Get sample name
  sample_name = fit$sample[1]
  
  # Fit plot
  colormap = ggsci::pal_jama("default")(7)[1:3]
  names(colormap) = c("Clonal", "Clonal LOH", "Subclonal")
  fit_plot = fit %>%
    ggplot2::ggplot() +
    CNAqc:::my_ggplot_theme() +
    ggplot2::geom_histogram(binwidth = 0.01, aes(VAF, fill = class)) +
    ggsci::scale_fill_jama() +
    xlim(0, 1) +
    ggplot2::labs(
      title = sample_name,
      subtitle = paste0("Purity: ", fit$purity[1], " - ", target_gene, ": ", gene_status)
    ) +
    ggplot2::guides(fill = guide_legend(''),
           color = guide_legend('', override.aes = aes(fill = NA)))  +
    ggplot2::geom_vline(aes(xintercept = gene_vafs,
                   color = "Target VAF"),
               linetype = 6) +
    ggplot2::geom_vline(
      aes(xintercept = fit$purity[1],
          color = "Purity"),
      linetype = 'dashed',
      size = .6
    ) +
    ggplot2::geom_vline(aes(xintercept = fit$purity[1] / 2,
                   color = "Clonal peak"),
               linetype = 'dashed') +
    ggplot2::scale_color_manual(
      name = "",
      values = c(
        `Target VAF` = "indianred3",
        `Purity` = "black",
        `Clonal Peak` = "gray"
      )
    )
  
  return(fit_plot)
}