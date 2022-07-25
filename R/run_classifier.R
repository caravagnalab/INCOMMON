#' Classify mutations using a (Beta-)Binomial model-based test.
#'
#' @param x An object of class \code{'TAPACLOTH'} generated with function `init`.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model Model used for the classification task, either "Binomial" (no over-dispersion) 
#' or "Beta-Binomial" (over-dispersion included), that will be used as the expected 
#' distribution for the number of reads with variant at fixed coverage and purity.
#' @param rho If "Beta-Binomial" model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#' @return An object of class `TAPACLOTH` containing the input plus
#' the classification data and paramters.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' x = init(mutations = example_data$data,
#'          sample = example_data$sample,
#'          purity = example_data$purity)
#' x = run_classifier(
#'     x, 
#'     alpha_level = 1e-3, 
#'    model = "Binomial")
#' print(x)
run_classifier = function(x,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
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
    
  pvalues = lapply(x$data$id, function(id) {
    null_model = test_setup(
      coverage = get_DP(x, mutation_id = id),
      purity = get_purity(x),
      rho = rho,
      alpha_level = alpha_level,
      model = model
    )
    
    pvalues = get_pvalues(x, null_model, id)
    pvalues$pvalue = p.adjust(pvalues$pvalue, method = "BH")
    pvalues$outcome = pvalues$pvalue < 1 - alpha_level
    
    return(pvalues)
  }) %>% 
    do.call(rbind, .)
  
  if ((model %>% tolower()) == "binomial") {
    test$classifier$binomial = list(
      params = tibble(alpha = alpha_level),
      data = full_join(x %>% idify() %>% get_data(), pvalues, by = c("id", "gene"))
    )
  }
  
  if ((model %>% tolower()) == "beta-binomial") {
    test$classifier$`beta-binomial` = list(
      params = tibble(alpha = alpha_level,
                      rho = rho),
      data = full_join(x %>% idify() %>% get_data(), pvalues, by = c("id", "gene"))
    )
  }
  return(test)
}


