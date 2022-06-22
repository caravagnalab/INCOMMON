plot_tapacloth = function(fit, target_gene){
  
  fit_plot = plot_fit(fit = fit$fit, target_gene = target_gene)
  
  target_data = fit$fit %>% 
    filter(gene == target_gene)
  
  target_plots = lapply(1:(target_data %>% nrow), function(i){
    
    null_model = test_setup(
    coverage = target_data$dp[i],
    purity = target_data$purity[i],
    rho = fit$rho,
    alpha_level = fit$alpha_level,
    model = fit$model
  )
    
  fit_power = plot_test_power(null_model)+
    geom_vline(xintercept = target_data$nv[i], linetype = 'dashed', size = .5)
  
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
