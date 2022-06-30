#' #' Plot for class \code{'TAPACLOTH'}.
#'

#' @description
#'
#' If the terzile model is used, the default plot is a histogram of the VAF spectrum is plotted,
#' with expected clonal peak, purity and VAF of mutations affecting the target gene
#' highlighted. If the Binomial or Beta-Binomial model is used, it includes
#' a plot showing details of the used statistical test in the bottom panel.
#'
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param target_gene The target gene for which classification of mutation calls is wanted.
#' @param sample_name The name of the sample to be analysed.
#' @param... Default S3 method parameter.
#'
#' @return Nothing.
#'
#' @export
#'
#' @examples
plot.TAPACLOTH = function(x, target_gene, sample_name, model, ...) {
  stopifnot(inherits(x, "TAPACLOTH"))

  # if(any(class(x$fit) == "bmix"))
  # {
  #   return(x$plot_bmix)
  # }

  fit_plot = plot_fit(x, target_gene, sample_name, model)

  if (model == "terzile") {
    return(fit_plot)
  }

  else{
    target_data = x$data$data %>%
      dplyr::filter(gene == target_gene)

    target_plots = lapply(1:(target_data %>% nrow), function(i) {
      null_model = test_setup(
        coverage = target_data$dp[i],
        purity = filter(x$data$purity, sample == sample_name)$purity,
        rho = x$classifier[[model]]$params$rho,
        alpha_level = x$classifier[[model]]$params$alpha,
        model = model
      )

      fit_power = plot_test_power(null_model) +
        ggplot2::geom_vline(xintercept = target_data$nv[i],
                            linetype = 'dashed',
                            size = .5)

    })

    # Fig assembly
    lp = append(list(fit_plot), target_plots)

    figure = ggpubr::ggarrange(plotlist = lp,
                               nrow = lp %>% length,
                               ncol = 1)

    return(figure)
  }
}

#' #' Print for class \code{'TAPACLOTH'}.
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

  if ("purity_estimate" %in% names(x))
  {
    for(s in names(x$purity_estimate)){
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Purity estimate using ',
          crayon::bgYellow(crayon::black("[ BMix ] ")),
          'with ',
          ifelse(s == 'bbinomial', "Beta-Binomial", "Binomial"),
          'model'
        )
      )
      print(full_join(x$purity_estimate[[s]]$purity, x$purity_estimate[[s]]$reliability, by = "sample") %>% as_tibble())
    }

    cli::cli_rule(paste('Purity estimate using ',
                        crayon::bgYellow(crayon::black("[ BMix ] "))))
    cat("\n")

    x$purity_estimate
  }
  #   cli::cli_rule(
  #     paste(
  #       crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
  #       'Purity estimate using BMix: '
  #     )
  #   )
  #   cat("\n")
  #   cli::cli_h3("BMix fit")
  #   cat("\n")
  #
  #   print(x$fit)
  #   cat("\n")
  #
  #   cli::cli_rule(
  #     paste(
  #       crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
  #       'Data: '
  #     )
  #   )
  #   print(x$data)
  # }
  # else
  # {
  if ("classifier" %in% names(x)) {

    for(s in names(x$classifier)){
      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ TAPACLOTH ] ")),
          'Test using {.field {s}} model',
          ifelse(s == 'bbinomial',
                 ' with overdispersion parameter {.field {x$classifier$bbinomial$params$rho}} and significance level {.field {x$classifier$bbinomial$params$alpha}}',
                 ifelse(s == 'binomial', 'with significance level {.field {x$classifier$bbinomial$params$alpha}}', '')),
          ''
        )
      )
      print(cbind(x$data$data, x$classifier[[s]]$data) %>% as_tibble())
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
