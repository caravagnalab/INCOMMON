#' Plot results of classification and statistical test.
#'
#' @param fit A fit object as output from `analyse_sample()` function.
#' @param target_gene The target gene for which classification of mutation calls is wanted.
#'
#' @return A figure assembly including (top) the VAF spectrum calls with expected clonal
#' peak, purity and VAF of mutations affecting the target gene are highlighted; and (bottom)
#' a plot showing details of the used statistical test.
#' @export
#'
#' @examples
plot_tapacloth = function(fit, target_gene){
  
  fit_plot = plot_fit(fit = fit$fit, target_gene = target_gene)
  
  target_data = fit$fit %>% 
    dplyr::filter(gene == target_gene)
  
  target_plots = lapply(1:(target_data %>% nrow), function(i){
    
    null_model = test_setup(
    coverage = target_data$dp[i],
    purity = target_data$purity[i],
    rho = fit$rho,
    alpha_level = fit$alpha_level,
    model = fit$model
  )
    
  fit_power = plot_test_power(null_model)+
    ggplot2::geom_vline(xintercept = target_data$nv[i], linetype = 'dashed', size = .5)
  
  }
  )
  
  # Fig assembly
  lp = append(list(fit_plot), target_plots)
  
  figure = ggpubr::ggarrange(
    plotlist = lp,
    nrow = lp %>% length,
    ncol = 1
  )
  
  return(figure)
}
