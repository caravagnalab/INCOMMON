#' Classify mutations using a (Beta-)Binomial model-based test.
#'
#' @param mutations a tibble with columns chromosome `chr`, start position `from`, end position `to`,
#'   reference `ref` and alternative `alt` alleles, coverage `DP`, number
#'   of reads with variant `NV`, variant allelic frequency `VAF` gene name `gene` 
#'   as Hugo Symbol.
#' @param sample sample name
#' @param purity sample purity.
#' The input can be prepared using function `init`.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model If set to "Binomial" or "Beta-Binomial": classification is run 
#' using a Binomial (no over-dispersion) or Beta-Binomial (over-dispersion
#' included) as expected distribution for the number of reads with variant at
#' fixed coverage and purity.
#' @param rho If Beta-Binomial model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#' @param tpanel A tibble assigning a role among oncogene, tumor suppressor gene (TSG)
#' and fusion to each gene, with columns `gene` for gene name and `gene_role` for the
#' role. Default is taken from COSMIC cancer gene census and stored in [data].
#'
#' @return An object of class `TAPACLOTH` that represents the classified input data.
#' @export
#'
#' @import dplyr
#'
#' @examples
#' x = run_classifier(
#'     example_data, 
#'     alpha_level = 1e-3, 
#'    model = "Binomial")
#' print(x)
run_classifier = function(mutations,
                          sample,
                          purity,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA,
                          tpanel = TAPACLOTH::cancer_gene_census)
{
  x = list(
    data = mutations,
    sample = sample,
    purity = purity
  )
  
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
  
  x$mutations = left_join(mutations, tpanel, by = "gene")
  
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
