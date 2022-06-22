#' Classify mutations in a sample.
#'
#' @param data A tibble containing mutations with sample name (sample), gene name
#' (gene), number of reads with variant (nv), coverage (dp), variant allele
#' frequency (VAF), and sample purity (purity) as columns.
#' @param sample_name The name of the sample to be analysed.
#' @param alpha_level The significance level to be used in hypothesis testing.
#' @param model The expected distribution for the number of reads with variant
#' at fixed coverage and purity, that can be chosen between Binomial (no over-dispersion)
#' or Beta-Binomial (over-dispersion included).
#' @param rho If Beta-Binomial model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the test.
#'
#' @return A list containing the input table reduced to the selected sample, with
#' the additional 'class' column accounting for the classification, the chosen
#' model, over-dispersion parameter (NA if the model is not Beta-Binomial) and
#' significance level used in the test.
#' @export
#'
#' @examples
analyse_sample = function(data,
                          sample_name,
                          alpha_level = 0.01,
                          model = "Binomial",
                          rho = NA)
  {

  cli::cli_h1(sample_name)

  sample_data = data %>%
    dplyr::filter(sample == sample_name)

  sample_data$class = sapply(1:(sample_data %>% nrow), function(i) {
    null_model = test_setup(
      coverage = sample_data$dp[i],
      purity = sample_data$purity[i],
      rho = rho,
      alpha_level = alpha_level,
      model = model
    )

    run_test(nv = sample_data$nv[i], null_model = null_model)

  })

  return(list(
    fit = sample_data,
    model = model,
    rho = rho,
    alpha_level = alpha_level
  ))

}
