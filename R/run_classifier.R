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
  
  if (inherits(x, "TAPACLOTH")) {
    test = x
    x = test$data
    if(!("classifier" %in% names(test))) test$classifier = list()
  }
  else{
    test$data = x
    test$classifier = list()
  }

  if ((model %>% tolower()) %in% c("binomial", "beta-binomial")) {
    
    x = lapply(unique(x$data$sample), function(s) {
      cli::cli_h1("TAPACLOTH {.field {model}} clonality/Zygosity testing for sample {.field {s}}")
      cat("\n")
      
      sample_data = x$data %>%
        dplyr::filter(sample == s)
      
      sample_purity = dplyr::filter(x$purity, sample == s)$purity
      
      cli::cli_alert_info("Computing null model distributions and p-values.")
      
      sample_data$cumprob = sapply(1:(sample_data %>% nrow), function(i) {
        null_model = test_setup(
          coverage = sample_data$dp[i],
          purity = sample_purity,
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
      
      sample_data$p_subclonal = p.adjust(sample_data$p_subclonal, method = "BH")    
      sample_data$p_loh = p.adjust(sample_data$p_loh, method = "BH")
      sample_data$class = case_when(
        sample_data$p_subclonal <= alpha_level ~ "Subclonal",
        sample_data$p_loh <= alpha_level ~ "Clonal LOH",
        TRUE ~ "Clonal")
      
      return(sample_data)
      
    }) %>%
      do.call(rbind, .)
    
    if ((model %>% tolower()) == "binomial") {
      test$classifier$binomial = list(
        params = tibble(alpha = alpha_level),
        data = x %>% select(class, starts_with("p_"))
      )
    }
    
    if ((model %>% tolower()) == "beta-binomial") {
      test$classifier$bbinomial = list(
        params = tibble(alpha = alpha_level,
                      rho = rho),
        data = x %>% select(class, starts_with("p_"))
      )
    }
  }

  else{
    # x = lapply(unique(x$data$gene), function(g) {
    #   cli::cli_h1(g)
    # 
    #   gene_data = x$data %>%
    #     dplyr::filter(gene == g)
    # 
    #   terziles = quantile(gene_data$VAF / gene_data$purity, probs = c(0, 1, 0.33)) %>%
    #     round(2)
    # 
    #   gene_data = gene_data %>%
    #     dplyr::mutate(
    #       class = case_when(
    #         VAF / x$puritypurity >= terziles[3] ~ "Clonal LOH",
    #         VAF / purity < terziles[3] ~ "Subclonal/Clonal"
    #       )
    #     )
    # }) %>%
    #   do.call(rbind, .)
    
    x = lapply(unique(x$data$sample), function(s) {
      cli::cli_h1(s)

      sample_data = x$data %>%
        dplyr::filter(sample == s)
      
      sample_purity = dplyr::filter(x$purity, sample == s)$purity

      terziles = quantile(sample_data$VAF / sample_purity, 
                          probs = c(0, 1, 0.33)) %>%
        round(2) %>% sort()

      sample_data = sample_data %>%
        dplyr::mutate(
          class = case_when(
            VAF / sample_purity >= terziles[3] ~ "Clonal LOH",
            VAF / sample_purity < terziles[3] ~ "Subclonal/Clonal"
          )
        )
    }) %>%
      do.call(rbind, .)
    
    test$classifier$terzile = list(
      data = x %>% select(class)
    )
  }

  # test$data = x
  # test$classifier = list(
  #   model = model,
  #   rho = rho,
  #   alpha_level = alpha_level
  # )
  

  return(test)

}
