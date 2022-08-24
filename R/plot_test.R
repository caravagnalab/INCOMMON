#' Visualize classification results for individual observations.
#'
#' @param x An object of class \code{'TAPACLOTH'} containing the classification results, as 
#' produced by function `run_classifier`.
#' @param id The id (format chr:from:to:ref:alt) of the mutation.
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
plot_test = function(x, id, model){
  stopifnot(inherits(x,"TAPACLOTH"))
  colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
  names_colors = expand.grid(names(colors), 1:2) %>% apply(1, paste, collapse = ' ')
  colors = sapply(names_colors, function(n)
    colors[strsplit(n, ' ')[[1]][1]])
  names(colors) = names_colors
  
  colors = c(colors, `out_of_sample` = 'gray')
  
  mdata = x %>% get_classifier(model) %>%
    get_data() %>%
    filter(id == !!id)
  
  dataset = pull(mdata, density)[[1]]
  
  dataset %>%
    ggplot() +
    geom_line(aes(x = NV, y = density), color = 'gainsboro',  size = .3) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    geom_point(data = dataset %>% maximise(),
               aes(x = NV,
                   y = density,
                   color = label),
               size = 1) +
    geom_line(data = dataset %>% maximise(),
              aes(x = NV,
                  y = density,
                  color = label),
              size = .3) +
    labs(
      title = paste0(
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
      "Alleles: ",
      mdata$ploidy,
      "N (Mutated: ",
      mdata$multiplicity,
      "N, Wild Type: ",
      mdata$wt,
      "N)\nPurity: ",
      get_purity(x),
      ' - cutoff: ',
      x %>% get_params(model = model) %>% pull(cutoff)
    ),
    caption = paste0("Classified using ", model, " model")) +
    # guides(color = 'none') +
    # geom_hline(data = dataset %>% group_by(label) %>% summarise(cutoff) %>% unique() %>% filter(label != "out of sample"),
    #     mapping = aes(yintercept = cutoff,
    #     color = label),
    #     linetype = 'dashed',
    #     size = .5
    #   )
    geom_vline(
      xintercept = get_NV(x, id),
      color = ifelse(is.na(mdata$ploidy), 'indianred3', 'forestgreen'),
      linetype = 'dashed',
      size = .5
    )
}