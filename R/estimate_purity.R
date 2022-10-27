#' Infer sample purity using a mixture of Binomial or Beta-Binomial distributions.
#'
#' @param x An object of class \code{'TAPACLOTH'} generated with function `init`.
#' @param eps In case data is fitted with a mixture of 2 distributions, sample purity
#' is estimated using only that whose peak is closer to the input purity estimate
#' than the value of this parameter. If both distributions satisfy this condition,
#' the one which is closest to the input purity estimate is used for the new estimate.
#'
#' @return An object of class `TAPACLOTH` containing the inferred purity.
#' @export
#'
#' @examples
#' input = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity)
#' out = estimate_purity(x = input, model = "binomial", eps = 0.01)
#' print(out)
estimate_purity = function(x,
                           eps = 0.01) 
  {
  
  stopifnot(inherits(x, "TAPACLOTH"))
  
  # Output
  test = x
  if (!("purity_estimate" %in% names(test)))
    test$purity_estimate = list()
  
  cli::cli_h1(
    "TAPACLOTH purity estimate of sample {.field {get_sample(x)}} using a Beta-binomial model"
  )
  cat("\n")
  
  # Prepare data for BMix
  input = data.frame(successes = get_data(x) %>% pull(NV),
                     trials = get_data(x) %>% pull(DP))
  # Return NA if NVs are less than 3
  if (nrow(input) <= 3) {
    cli::cli_alert("There are less than 4 SNVs: purity will not be estimated")
    return(test)
  }
  
  fit = BMix::bmixfit(input, K.Binomials = 0, K.BetaBinomials = 1:3)
  n_binomials = fit$K["BB"]
  peaks = sort(fit$BB.params["mu",]) %>% as.double()
  purity_bmix = purity_from_fit(
    n_binomials = n_binomials,
    peaks = peaks,
    purity = get_purity(x),
    eps = eps
  )
  plot_bmix = BMix::plot.bmix(fit, input)
  
  # Prepare output
  fit$data = input
  fit$plot_bmix = BMix::plot.bmix(fit, fit$data)
  fit$eps = eps
  fit$purity = min(purity_bmix, 1) %>% round(2)
  fit$reliability = ifelse(get_purity(x) == 0.0, NA, 1-sqrt(((get_purity(x)-fit$purity)/fit$purity)**2))
  
  test$purity_estimate = fit
  
  return(test)
}
