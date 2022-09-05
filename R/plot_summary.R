#'Plotting function for class \code{'TAPACLOTH'}.
#' @description
#' Generates a composite heatmap with summary statistics of classification over multiple samples.
#' @param x A list of obj of class \code{'TAPACLOTH'} containing a `classifier`.
#' @param model Model used for the classification, either "Binomial" or "Beta-Binomial".
#' @import pheatmap
#' @import CNAqc
#' @import ggplot2
#' @import ggsci 
#' @import ggridges 
#' @importFrom dplyr filter mutate rename select %>% 
#' @return An object of class \code{'pheatmap'}.
#' @export
#' @examples
#' x = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity)
#' x = run_classifier(x)
#' plot_heatmap(x)
plot_summary = function(x, model = 'binomial')
{
  stopifnot(inherits(x, "TAPACLOTH"))
  
  model = model %>% tolower()
  
  # Tibble 2 matrix
  t2m = function(x){
    m =  x[, 2:ncol(x)] %>% as.matrix()
    colnames(m) = names(x)[2:length(names(x))]
    rownames(m) = x$gene
    m
  }
  
  input_classes_matrix = x %>% 
    group_by(gene,label) %>% 
    summarise(n=n()) %>% 
    tidyr::pivot_wider(names_from = label, values_from = n, values_fill=0)
  
  matrix = input_classes_matrix %>% t2m()
  
  pheatmap::pheatmap(
    mat = matrix,
    main = "c",
    cluster_cols = FALSE,
    fontsize_row = 6,
    fontsize = 7,
    cellwidth = 20,
    cellheight = 10,
    filename = "./prova_pheatmap.pdf"
  ) 
  
  }