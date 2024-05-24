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
#' # First load example data
#' data(MSK_genomic_data)
#' data(MSK_clinical_data)
#' # Initialize the INCOMMON object (note the outputs to screen)
#' x = init(genomic_data = MSK_genomic_data, clinical_data = MSK_clinical_data)
#' # Run INCOMMON classification
#' x = classify(x = x, priors = pcawg_priors, entropy_cutoff = NULL, rho = 0.01)
#' # Use the genome interpreter to get intepreted classifications
#' x = genome_interpreter(x = x)
#' # Plot the fraction of class
#' x = genome_interpreter(x = x)
#' @importFrom dplyr filter mutate rename select %>%
#' @importFrom patchwork wrap_plots
plot_class_fraction = function(x, tumor_type, gene, ...){

  toplot = class_frequency(x, tumor_type, gene, ...)

  if('state' %in% colnames(toplot)){
    p = toplot %>%
      ggplot2::ggplot()+
      ggplot2::geom_bar(ggplot2::aes(x = '', y = frequency, fill = state), stat = 'identity')+
      scale_color_INCOMMON_class(aes = 'fill')+
      ggplot2::xlab('')+ggplot2::ylab('Fraction')+
      ggplot2::coord_flip()+
      CNAqc:::my_ggplot_theme(cex = .8)+
      ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON state'))+
      ggplot2::labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))
  }

  if('class' %in% colnames(toplot)){
    p2 = toplot %>%
      dplyr::mutate(class = dplyr::case_when(
        grepl('with LOH', class) ~ 'with LOH',
        grepl('without LOH', class) ~ 'without LOH',
        grepl('with AMP', class) ~ 'with AMP',
        grepl('without AMP', class) ~ 'without AMP',
        grepl('Tier-2', class) ~ 'Tier-2',
      )) %>%
      ggplot2::ggplot()+
      ggplot2::geom_bar(ggplot2::aes(x = '', y = frequency, fill = class), stat = 'identity')+
      scale_color_INCOMMON_high_level_class(aes = 'fill')+
      ggplot2::xlab('')+ggplot2::ylab('Fraction')+
      ggplot2::coord_flip()+
      CNAqc:::my_ggplot_theme(cex = .8)+
      ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON class'))+
      ggplot2::labs(title = paste0(gene, ' (', tumor_type, ')'), subtitle = paste0('N = ', unique(toplot$N)))

    p = patchwork::wrap_plots(p, p2, ncol = 1)
  }

  cli::cli_alert_info('The frequency of states {.field {toplot$state}} are {.field {round(toplot$frequency,2)}}')

  return(p)

}
