plot_test = function(x) {
  plotmodels = lapply(names(x$classifier), function(model) {
    plotlist = lapply(x$classifier[[model]]$data$gene %>% unique(), function(g) {
      
      ## Get model parameters
      alpha_level = x$classifier[[model]]$params$alpha
      rho = x$classifier[[model]]$params$rho
      ## Get NV of mutation under test
      nvtest = x$classifier[[model]]$data %>%
        filter(gene == g) %>%
        pull(NV) %>%
        unique()
      ## Get sample purity
      purity = x$purity
      ## Get data for tested gene and model used
      gdata = x$classifier[[model]]$data %>% dplyr::filter(gene == g)
      ## Format data for plot
      y = lapply(1:(gdata %>% nrow()),
                 function(i) {
                   dp = gdata[i,]$DP
                   k = gdata[i,]$karyotype
                   m = gdata[i,]$multiplicity
                   ploidy = stringr::str_split(k, pattern = ":")[[1]] %>% as.integer() %>% sum()
                   
                   if ((model %>% tolower()) == "binomial") {
                     p = dbinom(
                       x = 1:dp,
                       size = dp,
                       prob = m * purity / (2 * (1 - purity) + purity * ploidy)
                     )
                   }
                   if ((model %>% tolower()) == "beta-binomial") {
                     p = VGAM::dbetabinom(
                       x = 1:dp,
                       size = dp,
                       rho = rho,
                       prob = m * purity / (2 * (1 - purity) + purity * ploidy)
                     )
                   }
                   tibble(
                     nv = 1:dp,
                     p = p,
                     karyotype = gdata[i,]$karyotype,
                     multiplicity = gdata[i,]$multiplicity,
                     outcome = ifelse(gdata[i,]$pvalue > alpha_level, "PASS", "FAIL"),
                     l_a = gdata[i,]$l_a,
                     r_a = gdata[i,]$r_a
                   )
                 }) %>% do.call(rbind, .)
      y = y %>%
        mutate(class = ifelse(outcome == "FAIL", "FAIL", karyotype))
      ## Build plot
      plt = ggplot() +
        ggridges::geom_density_ridges(
          y %>% filter(nv > l_a & nv < r_a),
          mapping = aes(
            x = nv,
            y = karyotype,
            height = p,
            fill = class
          ),
          stat = "identity",
          alpha = 0.8
        )
      
      if ((model %>% tolower()) == "beta-binomial") {
        model_string = bquote("Test using Beta-Binomial model with " * rho * " = " *
                                .(rho))
      } else{
        model_string = bquote("Test using Binomial model")
      }
      
      plt + CNAqc:::my_ggplot_theme() +
        scale_fill_manual(
          values = c(
            "FAIL" = "#BEBEBE66",
            "1:0" = "steelblue",
            "1:1" = "#228B22CC",
            "2:0" =  "turquoise4",
            "2:1" = "#FFA500CC",
            "2:2" = "firebrick3"
          ),
          limits = c(unique(y$class)),
          guide = "none"
        ) +
        geom_vline(xintercept = nvtest, linetype = "longdash") +
        CNAqc:::my_ggplot_theme() +
        coord_cartesian(clip = "off") +
        ggplot2::labs(
          x = 'NV',
          y = "Pr(X = NV)",
          caption = model_string,
          title = paste0("Mutation", " ; Gene ", g),
          subtitle = bquote(
            "DP = " * .(unique(gdata$DP)) * '; ' * pi * ' = ' * .(purity) * ', ' ~ alpha *
              ' = ' * .(alpha_level)
          )
        )
    })
    names(plotlist) = x$classifier[[model]]$data$gene %>% unique()
    return(plotlist)
  })
  names(plotmodels) = names(x$classifier)
  for(model in names(x$classifier)) {
    x$classifier[[model]]$plot_test = plotmodels[[model]]
  }
  return(x)
}


