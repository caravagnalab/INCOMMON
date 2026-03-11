#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param k_max The maximum value of total copy number to be included in the model.
#' @param priors_k_m A dplyr::tibble or data frame with columns `gene`, `tumor_type`, `k`, `m`, `N` and `n`
#' and `p` indicating tumor-specific or pan-cancer (PANCA) prior probabilities.
#' @param priors_eta A dplyr::tibble or data frame with columns `tumor_type`,`mean_eta`, `var_eta`, `N`,`alpha_eta` and `beta_eta`
#' providing parameters of the Gamma prior distribution over the per copy sequencing rate.
#' @param purity_error The expected error on the input sample purity estimate.
#' @param num_cores The number of cores to use for parallel stan sampling.
#' @param iter_warmup The number of iterations of the stan warmup phase.
#' @param iter_sampling The number of iterations of the stan sampling  phase.
#' @param num_chains The number of MCMC chains to be run in parallel.
#' @return An object of class `INCOMMON` containing the original input plus
#' the classification data and parameters.
#' @export
#' @examples
#' # First load example data
#' data(MSK_genomic_data)
#' data(MSK_clinical_data)
#' data(priors_pcawg_hmf)
#' data(priors_eta)
#' # Initialize the INCOMMON object for a single sample (note the outputs to screen)
#' sample = 'P-0002081'
#' x = init(genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,], clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,])
#' # Run INCOMMON classification
#' x = classify(x = x, priors_k_m = priors_pcawg_hmf, priors_eta = priors_eta, num_cores = 1, iter_warmup = 10, iter_sampling = 10, num_chains = 1)
#' # An S3 method can be used to report to screen what is in the object
#' print(x)
#' @importFrom dplyr filter mutate rename select everything %>%

classify = function(
    x,
    k_max = 8,
    priors_k_m = priors_pcawg_hmf,
    priors_eta = priors_eta,
    purity_error = 0.05,
    num_cores = NULL,
    iter_warmup = 500,
    iter_sampling = 1000,
    num_chains = 4
    )
  {

  stopifnot(inherits(x, "INCOMMON"))

  if(length(samples(x))>1)
    cli::cli_abort("You can only fit one sample at a time, input with multiple samples provided!")

  check_input(x)

  cli::cli_alert_info("Performing classification")

  cli::cli_h1(
    "INCOMMON inference of copy number and mutation multiplicity for sample {.field {samples(x)}}"
  )
  cat("\n")

  ## Prepare input

  M = input(x) %>% nrow()
  N = input(x)$DP
  n = input(x)$NV
  purity_mean = purity(x = x, sample = samples(x))
  tumor_type = input(x) %>% dplyr::pull(tumor_type) %>% unique()

  priors_k_m_sample = get_stan_input_priors(x = x, N_mutations = M, priors = priors_k_m, k_max = k_max)
  priors_eta_sample = priors_eta %>% filter(tumor_type == !!tumor_type)
  if(nrow(priors_eta_sample)==0){
    priors_eta_sample = priors_eta %>% filter(tumor_type == 'PANCA')
  }

  data = list(
    M = M,
    N = N,
    n = n,
    k_max = k_max,
    purity_mean = purity_mean,
    purity_error = purity_error,
    alpha_x = priors_eta_sample$alpha_eta,
    beta_x = priors_eta_sample$beta_eta,
    alpha_k_m = priors_k_m_sample
  )

  ## Compile stan code

  model = get_stan_model()

  ## Fit model

  fit = model$sample(
    data = data,
    seed = 1992,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    chains = num_chains,
    parallel_chains = num_cores,
  )

  draws = draw_samples(fit, num_chains = num_chains, iter_sampling = iter_sampling, M = M, k_max = k_max)

  z_km = fit$summary(variables = 'z_km')

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% dplyr::arrange(k, m)

  z_km = lapply(1:M, function(i){
    z_km %>%
      dplyr::filter(grepl(paste0('z_km\\[',i,','), variable)) %>%
      dplyr::bind_cols(k_m_table) %>%
      dplyr::select(k, m, median) %>%
      dplyr::rename(z_km = median)
  })

  ## MAP estimates

  km_map = lapply(draws$km_rep, function(x){
    what = x %>%
      table() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(k = as.integer(k), m = as.integer(m))

    dplyr::arrange(what, dplyr::desc(n)) %>%
      dplyr::slice_head(n=1) %>%
      dplyr::select(-n)
  })

  eta_map = median(draws$eta_rep)
  purity_map = median(draws$purity_rep)

  ## Bayesian p-value tests

  bayes_p_eta = bayesian_p_value(posterior_rep = draws$eta_rep, prior_rep = draws$eta_prior_rep, prior_type = 'gamma')
  bayes_p_purity = bayesian_p_value(posterior_rep = draws$purity_rep, prior_rep = draws$purity_prior_rep, prior_type = 'beta')
  post_pred_DP = posterior_predictive_p_value(x = x, posterior_rep = draws$N_rep, observed_quantity = 'DP')
  post_pred_NV = posterior_predictive_p_value(x = x, posterior_rep = draws$n_rep, observed_quantity = 'NV')

  result = input(x) %>%
    dplyr::bind_cols(
      dplyr::tibble(
          purity_map,
          bayes_p_purity,
          eta_map,
          bayes_p_eta,
          post_pred_p.value_DP = post_pred_DP,
          post_pred_p.value_NV = post_pred_NV,
          z_km
        ) %>%
        dplyr::bind_cols(
          km_map %>% do.call(rbind, .) %>% dplyr::rename(map_k = k, map_m = m)
        )
    ) %>%
    dplyr::select(sample, tumor_type, purity, purity_map, eta_map, chr, from, to, gene, ref, alt, NV, DP, VAF, map_k, map_m, dplyr::everything())

  # if (!dir.exists(results_dir)) dir.create(results_dir, recursive = T)
  # if(!is.null(output)) saveRDS(object = output, file = paste0(results_dir,'/', sample, '.rds'))

  output = list(
    output = result,
    stan_fit = fit,
    parameters = dplyr::tibble(
      k_max,
      purity_error,
      stan_iter_warmup = iter_warmup,
      stan_iter_sampling = iter_sampling,
      num_chains
    ),
    priors_k_m = priors_k_m,
    priors_eta = priors_eta
    )

  class(output) = 'INCOMMON'

  output

}


