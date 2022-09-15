#' Visualize classification results for individual observations.
#'
#' @param x An object of class \code{'TAPACLOTH'} containing the classification results, as 
#' produced by function `run_classifier`.
#' @param model Model used for the classification task.
#' @return An object of class \code{'ggplot2'}.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity)
#' model = "beta-binomial"
#' x = run_classifier(
#'     x = x, 
#'     model = model, 
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' plot_test(x = x,id = x$classifier[[model]]data$id[1],model = model)
plot_test = function(x, model, assembly = F){
  stopifnot(inherits(x,"TAPACLOTH"))
  colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
  colors = colors[c("1:0","1:1","2:1","2:2")]
  names(colors) =sapply(names(colors), function(n){
    strsplit(n,split = ":")[[1]] %>% as.integer() %>% sum()
  })
  # names_colors = expand.grid(names(colors), 1:2) %>% apply(1, paste, collapse = ' ')
  # colors = sapply(names_colors, function(n)
  #   colors[strsplit(n, ' ')[[1]][1]])
  # names(colors) = sapply(names_colors, function(n){
  #   splitted = strsplit(n, ' ')[[1]]
  #   k = splitted[1]
  #   m = splitted[2]
  #   paste0(strsplit(k, ':')[[1]] %>% as.integer() %>% sum(),"N (Mutated: ", m,"N)")
  # })
  
  colors = c(colors, "out of sample" = 'gray')
  
  cutoff = x %>% get_params(model = model) %>% pull(cutoff)
  
  plots = lapply(x$classifier[[model]]$data$id, function(id){
    
    # uncertainty = round(x %>%  get_classifier() %>% get_data() %>% 
    #                       filter(id == !!id) %>% pull(uncertainty),2)
    
    
    mdata = x %>% get_classifier(model) %>%
      get_data() %>%
      filter(id == !!id)
    
    # dataset = pull(mdata, density)[[1]] %>%
    #   group_by(NV) %>%
    #   summarise(
    #     density,
    #     Major,
    #     minor,
    #     multiplicity,
    #     karyotype,
    #     label,
    #     peak,
    #     cutoff,
    #     uncertainty = 1 - max(density / sum(density))
    #   ) %>%
    #   ungroup()
    
    dataset = pull(mdata, density)[[1]]
    
    # label = dataset %>% 
    #   filter(NV==mdata$NV, ploidy == mdata$ploidy, multiplicity == mdata$multiplicity) %>% 
    #   pull(label)
    label = mdata$label
    
    scaleFactor = max(dataset$density)/max(dataset$uncertainty)
    
    dataset$multiplicity = factor(dataset$multiplicity)
    
    n_points = 1000
    stride = max(round(nrow(dataset)/n_points),1)
    
    dataset = dataset[seq(1,nrow(dataset),stride),] %>% 
      mutate(ploidy = ifelse(label=="out of sample", label, ploidy))
    
    dataset %>%
      ggplot() +
      geom_line(aes(x = NV, y = density), color = 'gainsboro',  size = .3)+
      geom_line(
        aes(x = NV, y = uncertainty * scaleFactor),
        color = "deeppink4",
        size = 0.5,
        linetype = "longdash",
        alpha = .4
      ) +
      scale_y_continuous("Likelihood", sec.axis = sec_axis(~./scaleFactor, name = "Uncertainty"))+
      CNAqc:::my_ggplot_theme() +
      theme(axis.title.y.right = element_text(color="deeppink4"),
            axis.text.y.right = element_text(color="deeppink4"))+
      # scale_color_manual(values = colors, breaks = dataset$label %>% unique()) +
      geom_point(data = dataset %>% maximise(),
                 aes(x = NV,
                     y = density,
                     color = ploidy,
                     shape = multiplicity
                     ),
                 size = 1) +
      scale_color_manual(values = colors)+
      scale_shape_manual(values = data_frame("1"=16,"2"=17))+
      geom_line(data = dataset %>% maximise(),
                aes(x = NV,
                    y = density,
                    color = ploidy,
                    size = multiplicity))+
      scale_size_manual(values = c("1"=.3,"2"=.3))+
      labs(
        title = paste0(
          mdata$from,
          "; ",
          mdata$ref,
          ">",
          mdata$alt,
          "; ",
          mdata$gene,
          " (",
          mdata$gene_role,
          ")"
        ),
        subtitle = paste0(
          "Classification: ",
          label,
          "\nSample: ",
          get_sample(x),
          "; Purity: ",
          get_purity(x)
        ),
        caption = bquote("Model: "*.(model)*";"~alpha*" = "*.(cutoff)~"; Uncertainty: "*.(round(mdata$uncertainty,2)))) +
      # guides(color = 'none') +
      # geom_hline(data = dataset %>% group_by(label) %>% summarise(cutoff) %>% unique() %>% filter(label != "out of sample"),
      #     mapping = aes(yintercept = cutoff,
      #     color = label),
      #     linetype = 'dashed',
      #     size = .5
      #   )
      geom_vline(
        xintercept = get_NV(x, id),
        # color = ifelse(is.na(mdata$ploidy), 'indianred3', 'forestgreen'),
        color = ifelse(label=="out of sample", 'indianred3', 'forestgreen'),
        linetype = 'dashed',
        size = .5
      )
  })
  if(assembly){
    ggpubr::ggarrange(plotlist = plots, ncol = 4)
  } else plots
}