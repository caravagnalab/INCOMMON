plot_fit = function(x, target_gene, sample_name, model) {
  # if(x$classifier$model == "terzile") {
  #   x$data = x$data %>%
  #     dplyr::mutate(class = class_terzile)
  # } else {
  #   x$data = x$data %>%
  #     dplyr::mutate(class = class_binom)
  # }

  # Select sample
  data = x$data$data %>%
    cbind(x$classifier[[ifelse((model %>% tolower())=="beta-binomial", "bbinomial", model %>% tolower())]]$data) %>%
    dplyr::filter(sample == sample_name)

  # Get class of mutations on target gene
  gene_status = data %>%
    dplyr::filter(gene == target_gene) %>%
    dplyr::mutate(class = paste0(class, ' (', nv, '/', dp, ')')) %>%
    dplyr::pull(class) %>%
    paste(collapse = ', ')

  # Get VAF of mutations on target gene
  gene_vafs = data %>%
    dplyr::filter(gene == target_gene) %>%
    dplyr::pull(VAF)

  # Fit plot
  colormap = ggsci::pal_jama("default")(7)[1:3]
  names(colormap) = c("Clonal", "Clonal LOH", "Subclonal")
  fit_plot = data %>%
    ggplot2::ggplot() +
    CNAqc:::my_ggplot_theme() +
    ggplot2::geom_histogram(binwidth = 0.01, aes(VAF, fill = class)) +
    ggplot2::scale_fill_manual(values = colormap)+
    # ggsci::scale_fill_jama() +
    xlim(0, 1) +
    ggplot2::labs(
      title = sample_name,
      subtitle = paste0("Purity: ", filter(x$data$purity, sample == sample_name)$purity, " - ", target_gene, ": ", gene_status)
    ) +
    ggplot2::guides(fill = guide_legend(''),
           color = guide_legend('', override.aes = aes(fill = NA)))  +
    ggplot2::geom_vline(aes(xintercept = gene_vafs,
                   color = "Target VAF"),
               linetype = 6) +
    ggplot2::geom_vline(
      aes(xintercept = filter(x$data$purity, sample == sample_name)$purity,
          color = "Purity"),
      linetype = 'dashed',
      size = .6
    ) +
    ggplot2::geom_vline(aes(xintercept = filter(x$data$purity, sample == sample_name)$purity / 2,
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