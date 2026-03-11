#' Generate report of the fit
#' @param x An object of class INCOMMON.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
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

