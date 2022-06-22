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
#' @return A list containing:
#' -the input table reduced to the selected sample, with
#' additional columns `purity_bmix` for the inferred purity and `purity_error` for
#' the relative error between the inferred and the input purity, assuming the 
#' inferred value is correct.
#' -an object of class `bmix` that represents a fit mixtur
#' -the inferred sample purity and 
#' -a cowplot figure showing results.
#' @export
#'
#' @examples
estimate_purity = function(data,
                           sample_name,
                           model = "Binomial",
                           purity = 1.0,
                           eps = 0.01) {

  cli::cli_h1(sample_name)

  sample_data = data %>%
    filter(sample == sample_name)
  
  # Prepare data for BMix
  nvs = sample_data$nv
  coverage = sample_data$dp
  input = data.frame(successes = nvs,
                     trials = coverage)

  # Return NA if NVs are less than 3
  if (length(nvs) <= 3)
    return(list(
      data = sample_data,
      fit = NA,
      purity = NA,
      plot_bmix = NA
    ))

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

  # Prepare output
  fit$data = input
  plot_bmix = BMix::plot.bmix(fit, fit$data)
  sample_data$purity_bmix = bmix_best_purity
  sample_data$purity_error = sqrt(((purity-bmix_best_purity)/bmix_best_purity)**2)

  return(list(
    data = sample_data,
    fit = fit,
    purity = bmix_best_purity,
    plot_bmix = plot_bmix
  ))
}
