#' Generate a diagnostic report of an INCOMMON model fit
#'
#' Combines, for a single sample, the posterior predictive checks for the
#' Poisson (depth) and Binomial (variant reads) sub-models, the eta and
#' purity posterior predictive checks, and the prior/posterior distributions
#' of copy number and multiplicity into one composite figure.
#'
#' @param x A classified object of class INCOMMON (see \code{\link{classify}}).
#' @return A composite \code{ggplot2}/\code{patchwork} object with the full
#' diagnostic report for the sample.
#' @export
#' @examples
#' \dontrun{
#' # plot_report requires a classified INCOMMON object (see ?classify)
#' data(MSK_genomic_data)
#' data(MSK_clinical_data)
#' data(priors_pcawg_hmf)
#' data(priors_eta)
#' sample = 'P-0002081'
#' x = init(
#'   genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,],
#'   clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,]
#' )
#' x = classify(
#'   x = x, priors_k_m = priors_pcawg_hmf, priors_eta = priors_eta,
#'   num_cores = 1, iter_warmup = 10, iter_sampling = 10, num_chains = 1
#' )
#' plot_report(x)
#' }
#' @importFrom dplyr filter mutate rename select %>% tibble intersect
#' @importFrom stats rgamma

plot_report = function(x){

  inp = input(x)

  num_chains = x$parameters$num_chains
  iter_sampling = x$parameters$stan_iter_sampling
  M = inp %>% nrow()
  k_max = x$parameters$k_max

  sample_id = samples(x)

  draws = draw_samples(
    fit = x$stan_fit,
    num_chains = x$parameters$num_chains,
    iter_sampling = x$parameters$stan_iter_sampling,
    M = x$output %>% nrow(),
    k_max = x$parameters$k_max
    )

  km_map = lapply(draws$km_rep, function(x){
    what = x %>%
      table() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(k = as.integer(k), m = as.integer(m))

    dplyr::arrange(what, dplyr::desc(n)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::select(-n)
  })

  p_pois = plot_poisson_model(
    x = x,
    sample = sample_id,
    N_rep = draws$N_rep,
    km_rep = draws$km_rep,
    km_map = km_map,
    purity_map = unique(x$output$purity_map),
    eta_map = unique(x$output$eta_map),
    post_pred_DP = x$output$post_pred_p.value_DP,
    k_max = k_max
    )

  p_eta = plot_eta_check(
    posterior_eta_rep = draws$eta_rep,
    prior_eta_rep = draws$eta_prior_rep,
    bayes_p = x$output$bayes_p_eta
    )

  p_purity = plot_purity_check(
    posterior_purity_rep = draws$purity_rep,
    prior_purity_rep = draws$purity_prior_rep,
    bayes_p = x$output$bayes_p_purity
    )

  p_eta_pi = patchwork::wrap_plots(p_eta, p_purity, design = 'A\nB', guides = 'collect')&
    ggplot2::theme(legend.position = 'bottom')

  pbin = plot_binomial_model(
    x = x,
    n_rep = draws$n_rep,
    km_rep = draws$km_rep,
    post_pred_NV = x$output$post_pred_p.value_NV
    )

  pkmprior = plot_prior_k_m(priors_k_m = x$priors_k_m, x = x, k_max = k_max)
  pkmpost = plot_posterior_k_m(x = x, k_max = k_max, z_km = x$output$z_km)

  plot_report = patchwork::wrap_plots(
    p_pois, p_eta_pi, pkmprior, pkmpost, pbin,
    design = 'AACDE\nBBCDE', guides = 'collect'
  ) +
    patchwork::plot_annotation(
      title = paste0(sample_id),
      subtitle = paste0(unique(x$output$tumor_type), '; M = ', M, ' mutations'))&ggplot2::theme(legend.position = 'bottom')

  plot_report
}

