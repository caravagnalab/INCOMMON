#' Plot for class \code{'TAPACLOTH'}.
#'

#' @description
#'
#' If the terzile model is used, the default plot is a histogram of the VAF spectrum is plotted,
#' with expected clonal peak, purity and VAF of mutations affecting the target gene
#' highlighted. If the Binomial or Beta-Binomial model is used, it includes
#' a plot showing details of the used statistical test in the bottom panel.
#'
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
# plot.TAPACLOTH = function(x, target_gene, sample_name, model, ...) {
#   stopifnot(inherits(x, "TAPACLOTH"))
# 
#   # if(any(class(x$fit) == "bmix"))
#   # {
#   #   return(x$plot_bmix)
#   # }
# 
#   fit_plot = plot_fit(x, target_gene, sample_name, model)
# 
#   if (model == "terzile") {
#     return(fit_plot)
#   }
# 
#   else{
#     target_data = x$data$data %>%
#       dplyr::filter(gene == target_gene)
# 
#     target_plots = lapply(1:(target_data %>% nrow), function(i) {
#       null_model = test_setup(
#         coverage = target_data$dp[i],
#         purity = filter(x$data$purity, sample == sample_name)$purity,
#         rho = x$classifier[[model]]$params$rho,
#         alpha_level = x$classifier[[model]]$params$alpha,
#         model = model
#       )
# 
#       fit_power = plot_test_power(null_model) +
#         ggplot2::geom_vline(xintercept = target_data$nv[i],
#                             linetype = 'dashed',
#                             size = .5)
# 
#     })
# 
#     # Fig assembly
#     lp = append(list(fit_plot), target_plots)
# 
#     figure = ggpubr::ggarrange(plotlist = lp,
#                                nrow = lp %>% length,
#                                ncol = 1)
# 
#     return(figure)
#   }
# }

plot.TAPACLOTH = function(x, ...) {
  stopifnot(inherits(x, "TAPACLOTH"))
  lapply(names(x$classifier), function(model){
    if("plot_test" %in% names(x$classifier[[model]])){
      return(x$classifier[[model]]$plot_test)
    }
    else{
      x = plot_test(x)
      return(x$classifier[[model]]$plot_test)
    }
  })
  
}
#' Print for class \code{'TAPACLOTH'}.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
print.TAPACLOTH = function(x, ...) {
  stopifnot(inherits(x, "TAPACLOTH"))
  ## Print input data
  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
      'Input data for sample {.field {x$sample}}, with purity {.field {x$purity}}'
    )
  )
  
  print(x$data)

  if ("purity_estimate" %in% names(x))
  {
    for(model in names(x$purity_estimate)){
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Purity estimate using ',
          crayon::bgYellow(crayon::black("[ BMix ] ")),
          'with {.field {model}} model'
        )
      )
      print(tibble(
        purity = get_purity(x),
        reliability = get_reliability(x, model),
        purity_bmix = get_purity_bmix(x, model),
      ))
    }
    
    cat("\n")
  }
  
  if ("classifier" %in% names(x)) {

    for(model in names(x$classifier)){
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Test using {.field {model}} model',
          ifelse((model %>% tolower()) == 'beta-binomial',
                 ' with overdispersion parameter {.field {x$classifier[[model]]$params$rho}} and significance level {.field {x$classifier[[model]]$params$alpha}}',
                 ' with significance level {.field {x$classifier[[model]]$params$alpha}}')),
          ''
        )
      print(get_classifier(x) %>% 
              get_data() %>% 
              dplyr::select(id, NV, DP, VAF, 
                            gene, gene_role, 
                            karyotype, multiplicity, 
                            l_a, r_a, pvalue, outcome)
            )
    }
  }
}

# Volevo fare qualcosa ma poi sono morto vedendo dentro x!
print_sample_test_binomial = function(x, sample)
{
  # x$fit %>% dplyr::filter()
  #
  # cli::cli_text("Classification: {.value {crayon:red()}}")

}
