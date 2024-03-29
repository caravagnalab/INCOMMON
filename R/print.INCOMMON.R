#' Print for class \code{'INCOMMON'}.
#' @param x An obj of class \code{'INCOMMON'}.
#' @param ... Default S3 method parameter.
#'
#' @return Nothing.
#' @import cli
#' @import crayon
#' @importFrom dplyr filter mutate rename select %>%
#' @export
print.INCOMMON = function(x, ...) {

  stopifnot(inherits(x, "INCOMMON"))

  n_pass_mutations = genomic_data(x, PASS = TRUE) %>% nrow()
  n_samples = genomic_data(x, PASS = TRUE) %>% dplyr::pull(sample) %>% unique() %>% length()
  n_genes = genomic_data(x, PASS = TRUE) %>% dplyr::pull(gene) %>% unique() %>% length()
  n_ttypes = clinical_data(x, PASS = TRUE) %>% dplyr::pull(tumor_type) %>% unique() %>% length()

  cli::cli_rule(
    paste(
      crayon::bgMagenta(crayon::black("[ INCOMMON ] ")),
      '{.field {n_pass_mutations}} PASS mutations across {.field {n_samples}} samples, with {.field {n_genes}} mutant genes across {.field {n_ttypes}} tumor types'
    )
  )


  mean_purity = clinical_data(x, PASS = TRUE) %>% dplyr::pull(purity) %>% mean(na.rm = T) %>% round(2)

  cli::cli_alert_info('Average sample purity: {.field {mean_purity}}')

  mean_depth = genomic_data(x, PASS = TRUE) %>% dplyr::pull(DP) %>% mean(na.rm = T) %>% as.integer()

  cli::cli_alert_info('Average sequencing depth: {.field {mean_depth}}')

  ## Print input data


  if ("classification" %in% names(x) & length(x$classification)!=0) {

      cli::cli_rule(
        paste(
          crayon::bgMagenta(crayon::black("[ INCOMMON ] ")),
          'Classified mutations with overdispersion parameter {.field {parameters(x)$rho}} and entropy cutoff {.field {parameters(x)$entropy_cutoff}}',
          ''
        ))
    
    
    print(classification(x))
  } else{

    print(x$input)
  }
}
