#' Classify mutations per sample using a (Beta-)Binomial model-based test, or
#' per gene using a 3-quantile test.
#'
#' @param x A list containing:
#' - `mutations`: a tibble with columns chromosome `chr`, start position `from`, end position `to`,
#'   reference `ref` and alternative `alt` alleles, coverage `DP`, number
#'   of reads with variant `NV`, variant allelic frequency `VAF` gene name `gene` as Hugo Symbol and
#'   `gene_role`.
#' - `sample` sample name
#' - `purity` sample purity.
#' The input can be prepared using function `init`.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model If set to "Binomial" or "Beta-Binomial": classification is run 
#' using a Binomial (no over-dispersion) or Beta-Binomial (over-dispersion
#' included) as expected distribution for the number of reads with variant at
#' fixed coverage and purity.
#' @param rho If Beta-Binomial model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#'
#' @return An object of class `TAPACLOTH` that represents the classified input data.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' data = list(data = dplyr::tibble(sample = "test", gene = c(paste("test gene ", 1:9), "target gene"), nv = c(seq(10,90,10), 120), dp = c(rep(100, 9), 200), VAF = c(seq(10,90,10), 120)/c(rep(100, 9), 200)), purity = dplyr::tibble(sample = "test", purity = 1))
#' data = run_classifier(x = data, alpha_level = 1e-3, model = "Binomial")
#' print(data)
run_classifier = function(x,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
{
  model = model %>% tolower()
  
  # Output
  test = list()
  class(test) = "TAPACLOTH"
  
  stopifnot(model%in% c("binomial", "beta-binomial", "terzile"))
  
  if (inherits(x, "TAPACLOTH")) {
    test = x
    if (!("classifier" %in% names(test)))
      test$classifier = list()
  }
  else{
    test = x
    test$classifier = list()
    class(test) = "TAPACLOTH"
  }
  
  if (model %in% c("binomial", "beta-binomial")) {
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
      pvalues$outcome = pvalues$pvalue > alpha_level
      
      return(pvalues)
    }) %>% do.call(rbind, .)
    
    if ((model %>% tolower()) == "binomial") {
      test$classifier$binomial = list(
        params = tibble(alpha = alpha_level),
        data = full_join(test %>% idify() %>% get_data(), pvalues, by = c("id", "gene"))
      )
    }
    
    if ((model %>% tolower()) == "beta-binomial") {
      test$classifier$`beta-binomial` = list(
        params = tibble(alpha = alpha_level,
                        rho = rho),
        data = full_join(test %>% idify() %>% get_data(), pvalues, by = c("id", "gene"))
      )
    }
  }
  
  # else{
  #   # x = lapply(unique(x$data$gene), function(g) {
  #   #   cli::cli_h1(g)
  #   #
  #   #   gene_data = x$data %>%
  #   #     dplyr::filter(gene == g)
  #   #
  #   #   terziles = quantile(gene_data$VAF / gene_data$purity, probs = c(0, 1, 0.33)) %>%
  #   #     round(2)
  #   #
  #   #   gene_data = gene_data %>%
  #   #     dplyr::mutate(
  #   #       class = case_when(
  #   #         VAF / x$puritypurity >= terziles[3] ~ "Clonal LOH",
  #   #         VAF / purity < terziles[3] ~ "Subclonal/Clonal"
  #   #       )
  #   #     )
  #   # }) %>%
  #   #   do.call(rbind, .)
  #
  #   x = lapply(unique(x$data$sample), function(s) {
  #     cli::cli_h1(s)
  #
  #     sample_data = x$data %>%
  #       dplyr::filter(sample == s)
  #
  #     sample_purity = dplyr::filter(x$purity, sample == s)$purity
  #
  #     terziles = quantile(sample_data$VAF / sample_purity,
  #                         probs = c(0, 1, 0.33)) %>%
  #       round(2) %>% sort()
  #
  #     sample_data = sample_data %>%
  #       dplyr::mutate(
  #         class = case_when(
  #           VAF / sample_purity >= terziles[3] ~ "Clonal LOH",
  #           VAF / sample_purity < terziles[3] ~ "Subclonal/Clonal"
  #         )
  #       )
  #   }) %>%
  #     do.call(rbind, .)
  #
  #   test$classifier$terzile = list(
  #     data = x %>% select(class)
  #   )
  # }
  
  # test$data = x
  # test$classifier = list(
  #   model = model,
  #   rho = rho,
  #   alpha_level = alpha_level
  # )
  
  
  return(test)
  
}
