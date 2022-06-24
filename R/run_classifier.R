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
run_classifier = function(data,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
{
  
  # Output
  test = list()
  class(test) = "TAPACLOTH"
  
  if (model %in% c("Binomial", "Beta-Binomial")) {
    data = lapply(unique(data$sample), function(s) {
      cli::cli_h1(s)
      
      sample_data = data %>%
        dplyr::filter(sample == s)
      
      sample_data$class_binom = sapply(1:(sample_data %>% nrow), function(i) {
        null_model = test_setup(
          coverage = sample_data$dp[i],
          purity = sample_data$purity[i],
          rho = rho,
          alpha_level = alpha_level,
          model = model
        )
        
        run_test(nv = sample_data$nv[i], null_model = null_model)
        
      })
      
      return(sample_data)
      
    }) %>%
      do.call(rbind, .)
  }
  
  else{
    data = lapply(unique(data$gene), function(g) {
      cli::cli_h1(g)
      
      gene_data = data %>%
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
  
  test$fit = data
  test$model = model
  test$rho = rho
  test$alpha_level = alpha_level
  
  return(test)
  
}
