#' Visualize classification results for individual observations.
#'
#' @param x An object of class \code{'INCOMMON'} containing the classification results, as 
#' produced by function `classify`.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' x = classify(
#'     x = x, 
#'     priors = NULL,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' plot_classification(x = x, id = ids(x)[1])
plot_classification = function(x, assembly = F){
  
  stopifnot(inherits(x,"INCOMMON"))
  
  plots = lapply(ids(x), function(id){
    
    # Get posterior of specific classified mutations in easy-to-plot format
    toplot = posterior(x, id = id)
    toplot = toplot %>% dplyr::group_by(NV) %>% dplyr::reframe(max = (value == max(value)), dplyr::across(dplyr::everything()))
      
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
      ggplot2::scale_size_manual(values = c("1"=.3,"2"=.3))+
      ggplot2::scale_y_continuous("Classification Probability", 
                         sec.axis = sec_axis(~./scaleFactor, name = "Classification Entropy"), 
                         limits = c(0,max(toplot$value)))+
      ggplot2::scale_color_manual(values = ploidy_colors)+
      CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.title.y.right = element_text(color="deeppink4"),
        axis.text.y.right = element_text(color="deeppink4"))+
      ggplot2::guides(size = "none",
             shape = ggplot2::guide_legend(title = 'Multiplicity', ncol = 1),
             color = ggplot2::guide_legend(title = 'Ploidy', ncol = 2))+
      ggplot2::labs(title = paste0(info(x, id)$gene, ' (', info(x, id)$state, ')'),
           subtitle = paste0(x$sample, ' (Purity = ', purity(x), '); ', info(x, id)$chr, '; ', info(x, id)$from, '; ', info(x, id)$ref, '>', info(x, id)$alt),
           caption = paste0('Entropy cut-off: ', parameters(x)$entropy_cutoff, '; Overdispersion: ', parameters(x)$rho))
      
  })
  if(assembly){
    ggpubr::ggarrange(plotlist = plots, ncol = 4)
  } else plots
}