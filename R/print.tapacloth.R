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
#' input = init(mutations = example_data$data[1,], sample = example_data$sample, purity = example_data$purity)
#' x = run_classifier(input, alpha_level = 1e-3, model = "Binomial")
#' x = plot_test(x)
#' plot.TAPACLOTH(x)

plot.TAPACLOTH = function(x, ...) {
  stopifnot(inherits(x, "TAPACLOTH"))
  lapply(names(x$classifier), function(model){
    lapply(x$classifier[[model]]$data$id, function(id){
      plot_test(x, id = id, model = model)
    })
  })
}
#' Print for class \code{'TAPACLOTH'}.
#' @param x An obj of class \code{'TAPACLOTH'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#' @import cli
#' @import crayon
#' @importFrom dplyr filter mutate rename select %>% 
#' @export
#'
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity)
#' x = run_classifier(
#'     x, 
#'     alpha_level = 1e-3, 
#'    model = "Binomial")
#' x = estimate_purity(x = x, model = "binomial", eps = 0.01)
#' print(x)
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
                 'with overdispersion parameter {.field {x$classifier[[model]]$params$rho}} and likelihood cutoff {.field {x$classifier[[model]]$params$cutoff}}',
                 'with likelihood cutoff {.field {x$classifier[[model]]$params$cutoff}}')),
          ''
        )
      print(get_classifier(x) %>% 
              get_data() %>% 
              dplyr::select(id, NV, DP, VAF, 
                            gene, gene_role, 
                            ploidy, multiplicity,wt,
                            uncertainty)
            )
    }
  }
}