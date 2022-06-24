#' Infer sample purity using a mixture of Binomial or Beta-Binomial distributions.
#'
#' @param data A tibble containing mutations with sample name (sample), gene name
#' (gene), number of reads with variant (nv), coverage (dp), variant allele
#' frequency (VAF), and sample purity (purity) as columns.
#' @param sample_name The name of the sample for which to infer purity.
#' @param model The expected distribution for the number of reads with variant
#' at fixed coverage and purity, that can be chosen between Binomial (no over-dispersion)
#' or Beta-Binomial (over-dispersion included).
#' @param purity An input purity estimate.
#' @param eps In case data is fitted with a mixture of 2 distributions, sample purity
#' is estimated using only that whose peak is closer to the input purity estimate
#' than the value of this parameter. If both distributions satisfy this condition,
#' the one which is closest to the input purity estimate is used for the new estimate.
#'
#' @return An object of class `TAPACLOTH` containing the inferred purity.
#' @export
#'
#' @examples
estimate_purity = function(data,
                           sample_name,
                           model = "Binomial",
                           purity = 1.0,
                           eps = 0.01) {
  # Output
  test = list()
  class(test) = "TAPACLOTH"
  
  cli::cli_h1(sample_name)

  sample_data = data %>%
    filter(sample == sample_name)
  
  if(is.na(purity)){
    cli::cli_alert("Input purity not available, reliability score will not be computed.")
    purity = 0.0
  }
  
  # Prepare data for BMix
  nvs = sample_data$nv
  coverage = sample_data$dp
  input = data.frame(successes = nvs,
                     trials = coverage)

  # Return NA if NVs are less than 3
  if (length(nvs) <= 3) {
    test$data = sample_data
    test$fit = NA
    test$purity = NA
    test$plot_bmix = NA
    return(test)
  }

  # Fit data with a mixture of 3 Binomials or BetaBinomials
  if (model == 'Binomial') {
    fit = BMix::bmixfit(input, K.Binomials = 1:3, K.BetaBinomials = 0)
    n_binomials = fit$K["B"]
    peaks = sort(fit$B.params)

  }

  else{
    fit = BMix::bmixfit(input, K.Binomials = 0, K.BetaBinomials = 1:3)
    n_binomials = fit$K["BB"]
    peaks = sort(fit$BB.params["mu", ]) %>% as.double()

  }
  # Estimate purity from fits when feasible
  bmix_best_purity = case_when(
    n_binomials == 3 ~ peaks[2] * 2,
    n_binomials == 1 ~ peaks[1] * 2,
    n_binomials == 2 &
      min(abs(peaks[2] * 2 - purity), abs(peaks[1] * 2 - purity)) <= eps ~ peaks[which.min(c(abs(peaks[2] *
                                                                                                   2 - purity), abs(peaks[1] * 2 - purity)))] * 2,
    TRUE ~ NA %>% as.double()

  )
  
  if(is.na(bmix_best_purity)){
    cli::cli_alert("Purity could not be estimated reliably. Check results of the fit in the output figure")
  }

  # Prepare output
  fit$data = input
  plot_bmix = BMix::plot.bmix(fit, fit$data)
  sample_data$purity_bmix = min(bmix_best_purity, 1) %>% round(2)
  sample_data$reliability = ifelse(purity == 0.0, NA, 1-sqrt(((purity-bmix_best_purity)/bmix_best_purity)**2))

  test$data = sample_data
  test$fit = fit
  test$purity = bmix_best_purity
  test$plot_bmix = plot_bmix
  
  return(test)
}
