#' Classify mutations using a (Beta-)Binomial model-based test.
#'
#' @param x An object of class \code{'TAPACLOTH'} generated with function `init`.
#' @param cutoff Likelihood cut-off for class assignment.
#' @param model Model used for the classification task, either "Binomial" (no over-dispersion) 
#' or "Beta-Binomial" (over-dispersion included), that will be used as the expected 
#' distribution for the number of reads with variant at fixed coverage and purity.
#' @param rho If "Beta-Binomial" model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @return An object of class `TAPACLOTH` containing the input plus
#' the classification data and paramters.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity)
#' x = run_classifier(
#'     x = x, 
#'     model = "Beta-Binomial", 
#'     rho = 0.01,
#'     karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#'     )
#' print(x)
run_classifier = function(x,
                          cutoff = 0.75,
                          model = "Binomial",
                          rho = 0.01,
                          karyotypes = c("1:0","1:1","2:0","2:1","2:2")
                          )
  {
  stopifnot(inherits(x, "TAPACLOTH"))
  
  model = model %>% tolower()
  stopifnot(model%in% c("binomial", "beta-binomial"))
  
  # Output
  test = x
  if (!("classifier" %in% names(test)))
    test$classifier = list()
  
  cli::cli_h1(
      "TAPACLOTH {.field {model}} clonality/Zygosity testing for sample {.field {x$sample}}"
    )
    cat("\n")
    
  cli::cli_alert_info("Computing null model distributions and p-values.")
    
  x = idify(x)
    
  tests = lapply(x$data$id, function(id) {
    
    binomial_test(
      test = get_NV(x, id),
      DP = get_DP(x, id),
      purity = get_purity(x),
      cutoff = cutoff,
      model = model,
      rho = rho,
      karyotypes = karyotypes,
    )
  }) %>% 
    do.call(rbind, .)
  
  if (model == "binomial") {
    test$classifier$binomial = list(
      params = tibble(cutoff = cutoff),
      data = bind_cols(x %>% get_data(), tests)
    )
  }
  
  if (model == "beta-binomial") {
    test$classifier$`beta-binomial` = list(
      params = tibble(cutoff = cutoff,
                      rho = rho),
      data = bind_cols(x %>% get_data(), tests)
    )
  }
  return(test)
}


