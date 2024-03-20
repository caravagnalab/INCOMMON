#' Print for class \code{'INCOMMON'}.
#' @param x An obj of class \code{'INCOMMON'}.
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
#'          purity = example_data$purity,
#'          tumor_type = example_data$tumor_type)
#' print(x)
#' x = classify(
#'     x = x, 
#'     priors = NULL,
#'     entropy_cutoff = 0.2,
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' print(x)
print.INCOMMON = function(x, ...) {
  stopifnot(inherits(x, "INCOMMON"))
  ## Print input data
  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::black("[ INCOMMON ] ")),
      'Sample {.field {x$sample}} ({.field {x$tumor_type}}), with purity {.field {x$purity}}'
    )
  )
  
  if ("fit" %in% names(x)) {

      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ INCOMMON ] ")),
          'Classified mutations using Beta-binomial model',
                 'with overdispersion parameter {.field {parameters(x)$rho}} and entropy cutoff {.field {parameters(x)$entropy_cutoff}}',
          ''
        ))
      print(classification(x))
  } else{
    print(data(x))
  }
}
