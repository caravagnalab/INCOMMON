#' Classify mutations per sample using a (Beta-)Binomial model-based test, or
#' per gene using a 3-quantile test.
#'
#' @param data A tibble containing mutations with sample name (sample), gene name
#' (gene), number of reads with variant (nv), coverage (dp), variant allele
#' frequency (VAF), and sample purity (purity) as columns.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model If set to "Binomial" or "Beta-Binomial": classification is run per
#' sample, using a Binomial (no over-dispersion) or Beta-Binomial (over-dispersion
#' included) as expected distribution for the number of reads with variant at
#' fixed coverage and purity. If set to "terzile": classification is run per gene,
#' and is based on a terzile test.
#' @param rho If Beta-Binomial model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#'
#' @return An object of class `TAPACLOTH` that represents the classified input data.
#' @export
#'
#' @import dplyr
#'
#' @examples
run_classifier = function(x,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
{
  # Output
  test = list()
  class(test) = "TAPACLOTH"

  stopifnot((model %>% tolower()) %in% c("binomial", "beta-binomial", "terzile"))

  if (model %in% c("Binomial", "Beta-Binomial")) {
    
    x = lapply(unique(x$sample), function(s) {
      cli::cli_h1("TAPACLOTH {.field {model}} clonality/Zygosity testing for sample {.field {s}}")
      cat("\n")
      
      sample_data = x %>%
        dplyr::filter(sample == s)
      
      cli::cli_alert_info("Computing null model distributions and p-values.")
      
      sample_data$cumprob = sapply(1:(sample_data %>% nrow), function(i) {
        null_model = test_setup(
          coverage = sample_data$dp[i],
          purity = sample_data$purity[i],
          rho = rho,
          alpha_level = alpha_level,
          model = model
        )
        
        # class = run_test(nv = sample_data$nv[i], null_model = null_model)
        cumprob = null_model$density$p[sample_data$nv[i]]
        return(cumprob)
      })
      
      sample_data = sample_data %>%
        dplyr::mutate(
          p_subclonal = cumprob,
          p_loh = 1 - cumprob,
          ) %>%
        select(-cumprob)
      
      
      return(sample_data)
      
    }) %>%
      do.call(rbind, .)
    
    x$p_subclonal = p.adjust(x$p_subclonal, method = "BH")    
    x$p_loh = p.adjust(x$p_loh, method = "BH")
    x$class_binom = case_when(
      x$p_subclonal <= alpha_level ~ "Subclonal",
      x$p_loh <= alpha_level ~ "Clonal LOH",
      TRUE ~ "Clonal")
      
  }

  else{
    x = lapply(unique(x$gene), function(g) {
      cli::cli_h1(g)

      gene_data = x %>%
        dplyr::filter(gene == g)

      terziles = quantile(gene_data$VAF / gene_data$purity, probs = c(0, 1, 0.33)) %>%
        round(2)

      gene_data = gene_data %>%
        dplyr::mutate(
          class_terzile = case_when(
            VAF / purity >= terziles[3] ~ "Clonal LOH",
            VAF / purity < terziles[3] ~ "Subclonal/Clonal"
          )
        )
    }) %>%
      do.call(rbind, .)
  }

  test$fit = x
  test$model = model
  test$rho = rho
  test$alpha_level = alpha_level

  return(test)

}
