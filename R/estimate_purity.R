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
estimate_purity = function(x,
                           model = "Binomial",
                           eps = 0.01) {
  
  # Output
  test = list()
  class(test) = "TAPACLOTH"
  
  stopifnot((model %>% tolower()) %in% c("binomial", "beta-binomial"))
  
  if (inherits(x, "TAPACLOTH")) {
    test = x
    x = test$data
    if(!("purity_estimate" %in% names(test))) test$purity_estimate = list()
  }
  else{
    test$data = x
    test$purity_estimate = list()
  }
  
  samples = unique(x$data$sample)
  
  x = lapply(samples, function(s) {
    cli::cli_h1("TAPACLOTH {.field {model}} clonality/Zygosity testing for sample {.field {s}}")
    cat("\n")
    
    sample_data = x$data %>%
      dplyr::filter(sample == s)
    
    sample_purity = dplyr::filter(x$purity, sample == s)$purity
    
    cli::cli_h1("TAPACLOTH purity estimate of sample {.field {s}} using {.field {model}} model")
    cat("\n")
  
  if(is.na(sample_purity)){
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
    cli::cli_alert("There are less than 3 SNVs: purity will not be estimated")
    return(test)
  }

  # Fit data with a mixture of 3 Binomials or BetaBinomials
  if (model == 'Binomial') {
    fit = BMix::bmixfit(input, K.Binomials = 1:3, K.BetaBinomials = 0)
    n_binomials = fit$K["B"]
    peaks = sort(fit$B.params)
    purity_bmix = purity_from_fit(n_binomials, peaks, sample_purity, eps)
  }

  else{
    fit = BMix::bmixfit(input, K.Binomials = 0, K.BetaBinomials = 1:3)
    n_binomials = fit$K["BB"]
    peaks = sort(fit$BB.params["mu", ]) %>% as.double()
    purity_bmix = purity_from_fit(n_binomials, peaks, sample_purity, eps)
    plot_bmix = BMix::plot.bmix(fit, input)
    # test$purity_estimate$binomial = list(params = tibble(alpha = alpha_level))
  }
  
  # Prepare output
  fit$data = input
  fit$plot_bmix = BMix::plot.bmix(fit, fit$data)
  fit$eps = eps
  fit$purity = min(purity_bmix, 1) %>% round(2)
  fit$reliability = ifelse(sample_purity == 0.0, NA, 1-sqrt(((sample_purity-fit$purity)/fit$purity)**2))
  return(fit)
  })
  
  plot_bmix = lapply(1:length(x), function(n){x[[n]]$plot_bmix})
  names(plot_bmix) = samples
  
    test$purity_estimate[[ifelse(model == "binomial", "binomial", "bbinomial")]] = list(
      params = tibble(eps = eps),
      purity = tibble(sample = samples,
                      purity = sapply(1:length(x), function(n){x[[n]]$purity})),
      reliability = tibble(sample = samples,
                           reliability = sapply(1:length(x), function(n){x[[n]]$reliability})),
      plot_bmix = plot_bmix
      )
    
  # test$purity_estimate = lapply(1:length(x), function(n){x[[n]]$fit})
  # names(test$purity_estimate) = samples
  
  # test$data = sample_data
  # test$fit = fit
  # test$purity = bmix_best_purity
  # test$plot_bmix = plot_bmix
  
  return(test)
}
