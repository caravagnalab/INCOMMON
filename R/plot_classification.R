#' Visualize classification results for individual observations.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results, as
#' produced by function `classify`.
#' @param sample Sample name from the dataset.
#' @param assembly Whether to assemble plots of mutations from the same sample in
#' a multi-faceted plot.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Plot classification results for a specific sample
#' plot_classification(x = MSK_classified, sample = 'P-0002081-T01-IM3')
#' @importFrom dplyr filter mutate rename select %>%
plot_classification = function(x, sample, assembly = F){

  stopifnot(inherits(x,"INCOMMON"))

  x = subset_sample(x, sample = sample)

  plots = lapply(ids(x), function(id){

    # Get posterior of specific classified mutations in easy-to-plot format
    toplot = posterior(x, id = id)
    toplot = toplot %>%
      dplyr::group_by(NV) %>%
      dplyr::reframe(max = (value == max(value)), dplyr::across(dplyr::everything()))

    toplot$multiplicity = factor(toplot$multiplicity)
    toplot$ploidy = factor(toplot$ploidy)

    highlight = toplot %>% dplyr::filter(max)
    highlight = lapply(split(highlight, highlight$label), function(s){
      if(nrow(s) >= 20){
        stride = as.integer(nrow(s)/20)
        s = s[seq(1,nrow(s),stride),]
      }
      s
    }) %>% do.call(rbind, .)

    # Get scale factor for overlapped visualisation of posterior and entropy
    scaleFactor = max(toplot$value)/max(toplot$entropy)
    scaleFactor = ifelse(is.infinite(scaleFactor), 1, scaleFactor)

    plot = toplot %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = NV, y = value),
                color = 'white',
                size = .5,
                alpha = 0)+
      ggplot2::geom_line(
        ggplot2::aes(x = NV, y = entropy * scaleFactor),
        color = "deeppink4",
        size = 0.5,
        linetype = "longdash",
        alpha = .4
      ) +
      ggplot2::geom_point(data = highlight,
                          ggplot2::aes(x = NV,
                     y = value,
                     color = ploidy,
                     shape = multiplicity
                 ),
                 size = 1) +
      ggplot2::geom_line(data = highlight,
                         ggplot2::aes(x = NV,
                     y = value,
                     color = ploidy,
                     group = label
                 ),
                 size = .5) +
      ggplot2::geom_vline(
        xintercept = NV(x, id),
        color = 'black',
        linetype = 'dashed',
        size = .5
      )+
      ggplot2::geom_hline(
        yintercept = entropy(x, id)*scaleFactor,
        color = 'deeppink4',
        linetype = 'dashed',
        size = .5
      )+
      ggplot2::scale_size_manual(values = c("1"=.3,"2"=.3))+
      ggplot2::scale_y_continuous("Classification Probability",
                         sec.axis = ggplot2::sec_axis(~./scaleFactor, name = "Classification Entropy"),
                         limits = c(0,max(toplot$value)))+
      ggplot2::scale_color_manual(values = ploidy_colors)+
      CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.title.y.right = ggplot2::element_text(color="deeppink4"),
        axis.text.y.right = ggplot2::element_text(color="deeppink4"))+
      ggplot2::guides(size = "none",
             shape = ggplot2::guide_legend(title = 'Multiplicity', ncol = 1),
             color = ggplot2::guide_legend(title = 'Ploidy', ncol = 2))+
      ggplot2::labs(
        title = paste0(info(x, id)$gene, ' (', info(x, id)$state, ')'),
        subtitle = paste0(x$sample, ' (Purity = ', purity(x, id), '); ', info(x, id)$chr, '; ', info(x, id)$from, '; ', info(x, id)$ref, '>', info(x, id)$alt),
        caption = paste0('Entropy cut-off: ', parameters(x)$entropy_cutoff, '; Overdispersion: ', parameters(x)$rho))

  })
  if(assembly){

    toplot = lapply(ids(x), function(id){posterior(x, id = id) %>% dplyr::mutate(gene = gene(x, id))}) %>%
                      do.call(rbind, .)
    toplot = toplot %>%
      dplyr::group_by(gene, NV) %>%
      dplyr::reframe(max = (value == max(value)), dplyr::across(dplyr::everything()))

    toplot$multiplicity = factor(toplot$multiplicity)
    toplot$ploidy = factor(toplot$ploidy)

    highlight = toplot %>% dplyr::filter(max)
    highlight = lapply(split(highlight, highlight$label), function(s){
      if(nrow(s) >= 20){
        stride = as.integer(nrow(s)/20)
        s = s[seq(1,nrow(s),stride),]
      }
      s
    }) %>% do.call(rbind, .)

    # Get scale factor for overlapped visualisation of posterior and entropy
    scaleFactor = max(toplot$value)/max(toplot$entropy)
    scaleFactor = ifelse(is.infinite(scaleFactor), 1, scaleFactor)

    plot = toplot %>%
      ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = NV, y = value),
                         color = 'white',
                         size = .5,
                         alpha = 0)+
      ggplot2::geom_line(
        ggplot2::aes(x = NV, y = entropy * scaleFactor),
        color = "deeppink4",
        size = 0.5,
        linetype = "longdash",
        alpha = .4
      ) +
      ggplot2::geom_point(data = highlight,
                          ggplot2::aes(x = NV,
                                       y = value,
                                       color = ploidy,
                                       shape = multiplicity
                          ),
                          size = 1) +
      ggplot2::geom_line(data = highlight,
                         ggplot2::aes(x = NV,
                                      y = value,
                                      color = ploidy,
                                      group = label
                         ),
                         size = .5) +
      ggplot2::geom_vline(
        data = lapply(ids(x), function(id) dplyr::tibble(id = id, gene = gene(x, id), NV = NV(x, id))) %>% do.call(rbind, .),
        ggplot2::aes(xintercept = NV),
        color = 'black',
        linetype = 'dashed',
        size = .5
      )+
      ggplot2::geom_hline(
        data = lapply(ids(x), function(id) dplyr::tibble(id = id, gene = gene(x, id), entropy = entropy(x, id)*scaleFactor)) %>% do.call(rbind, .),
        ggplot2::aes(yintercept = entropy),
        color = 'deeppink4',
        linetype = 'dashed',
        size = .5
      )+
      ggplot2::scale_size_manual(values = c("1"=.3,"2"=.3))+
      ggplot2::scale_y_continuous("Classification Probability",
                                  sec.axis = ggplot2::sec_axis(~./scaleFactor, name = "Classification Entropy"),
                                  limits = c(0,max(toplot$value)))+
      ggplot2::scale_color_manual(values = ploidy_colors)+
      ggplot2::facet_wrap(~gene) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.title.y.right = ggplot2::element_text(color="deeppink4"),
        axis.text.y.right = ggplot2::element_text(color="deeppink4"))+
      ggplot2::guides(size = "none",
                      shape = ggplot2::guide_legend(title = 'Multiplicity', ncol = 1),
                      color = ggplot2::guide_legend(title = 'Ploidy', ncol = 2))
    return(plot)
  } else plots
}
