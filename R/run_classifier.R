#' Classify mutations. The classification task can be run either per sample,
#' using either a Binomial or Beta-Binomial model-based statistical test, or
#' per gene, in which case the classification is based on a terzile test across
#' all samples.
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
#' @return A list containing the input table reduced to the selected sample, with
#' the additional 'class' column accounting for the classification, the chosen
#' model, over-dispersion parameter (NA if the model is not Beta-Binomial) and
#' significance level used in the test.
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
        round(1)
      
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
  
  return(list(
    fit = data,
    model = model,
    rho = rho,
    alpha_level = alpha_level
  ))
  
}
