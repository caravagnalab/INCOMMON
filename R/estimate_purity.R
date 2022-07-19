#' Infer sample purity using a mixture of Binomial or Beta-Binomial distributions.
#'
#' @param mutations a tibble with columns chromosome `chr`, start position `from`, end position `to`,
#'   reference `ref` and alternative `alt` alleles, coverage `DP`, number
#'   of reads with variant `NV`, variant allelic frequency `VAF` gene name `gene` 
#'   as Hugo Symbol.
#' @param sample sample name
#' @param purity sample purity.
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
#' data = list(data = dplyr::tibble(sample = "test", gene = paste0("test gene ", 1:30), nv =  c(seq(5, 14, 1), seq(40,58,2), seq(80, 98, 2))*2, dp = 200, VAF = c(seq(5, 14, 1), seq(40,58,2), seq(80, 98, 2))*2/200), purity = dplyr::tibble(sample = "test",  purity = 0.4))
#' data = estimate_purity(x = data, model = "binomial", eps = 0.01)
#' print(data)
estimate_purity = function(x,
                           model = "Binomial",
                           eps = 0.01,
                           tpanel = TAPACLOTH::cancer_gene_census) 
  {
  
  stopifnot(inherits(x, "TAPACLOTH"))
  
  model = model %>% tolower()
  stopifnot(model %in% c("binomial", "beta-binomial"))
  
  # Output
  test = x
  if (!("purity_estimate" %in% names(test)))
    test$purity_estimate = list()
  
  cli::cli_h1(
    "TAPACLOTH purity estimate of sample {.field {get_sample(x)}} using {.field {model}} model"
  )
  cat("\n")
  
  # Prepare data for BMix
  input = data.frame(successes = get_data(x) %>% pull(NV),
                     trials = get_data(x) %>% pull(DP))
  # Return NA if NVs are less than 3
  if (nrow(input) <= 3) {
    cli::cli_alert("There are less than 3 SNVs: purity will not be estimated")
    return(test)
  }
  # Fit data with a mixture of 3 Binomials or BetaBinomials
  if (model == 'binomial') {
    fit = BMix::bmixfit(input, K.Binomials = 1:3, K.BetaBinomials = 0)
    n_binomials = fit$K["B"]
    peaks = sort(fit$B.params)
    purity_bmix = purity_from_fit(
      n_binomials = n_binomials,
      peaks = peaks,
      purity = get_purity(x),
      eps = eps
    )
  }
  else{
    fit = BMix::bmixfit(input, K.Binomials = 0, K.BetaBinomials = 1:3)
    n_binomials = fit$K["BB"]
    peaks = sort(fit$BB.params["mu", ]) %>% as.double()
    purity_bmix = purity_from_fit(n_binomials = n_binomials, 
                                  peaks = peaks, 
                                  purity = get_purity(x), 
                                  eps = eps)
    plot_bmix = BMix::plot.bmix(fit, input)
    # test$purity_estimate$binomial = list(params = tibble(alpha = alpha_level))
  }
  
  # Prepare output
  fit$data = input
  fit$plot_bmix = BMix::plot.bmix(fit, fit$data)
  fit$eps = eps
  fit$purity = min(purity_bmix, 1) %>% round(2)
  fit$reliability = ifelse(get_purity(x) == 0.0, NA, 1-sqrt(((get_purity(x)-fit$purity)/fit$purity)**2))
  
  if ((model %>% tolower()) == "binomial") {
    test$purity_estimate$`binomial` = fit
  }
  else{
    test$purity_estimate$`beta-binomial` = fit
  }
  return(test)

  # plot_bmix = lapply(1:length(x), function(n){x[[n]]$plot_bmix})
  # names(plot_bmix) = samples
  # 
  #   test$purity_estimate[[model]] = list(
  #     params = tibble(eps = eps),
  #     purity = tibble(sample = samples,
  #                     purity = sapply(1:length(x), function(n){x[[n]]$purity})),
  #     reliability = tibble(sample = samples,
  #                          reliability = sapply(1:length(x), function(n){x[[n]]$reliability})),
  #     plot_bmix = plot_bmix
  #     )
    
  # test$purity_estimate = lapply(1:length(x), function(n){x[[n]]$fit})
  # names(test$purity_estimate) = samples
  
  # test$data = sample_data
  # test$fit = fit
  # test$purity = bmix_best_purity
  # test$plot_bmix = plot_bmix
  
  # return(test)
}
