#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param k_max The maximum value of total copy number to be included in the model.
#' @param priors_m_k A dplyr::tibble or data frame with columns `gene`, `gene_role`, `tumor_type`, `ploidy`, `multiplicity`
#' and `p` indicating tumor-specific or pan-cancer (PANCA) prior probabilities.
#' @param priors_x A dplyr::tibble or data frame with columns `mean_x`, `sigma2x`, `alpha_x` and `beta_x`
#' providing parameters of the Gamma prior distribution over the per copy sequencing rate.
#' @param purity_error The expected error on the input sample purity estimate.
#' @param num_cores The number of cores to use for parallel stan sampling.
#' @param  iter_warmup The number of iterations of the stan warmup phase.
#' @param iter_sampling The number of iterations of the stan sampling  phase.
#' @param dump Whether to dump results to a file.
#' @param dump_file The file path for dumping results.
#' @return An object of class `INCOMMON` containing the original input plus
#' the classification data and parameters.
#' @export
#' @examples
#' # First load example data
#' data(MSK_genomic_data)
#' data(MSK_clinical_data)
#' # Initialize the INCOMMON object for a single sample (note the outputs to screen)
#' sample = 'P-0002081'
#' x = init(genomic_data = MSK_genomic_data[MSK_genomic_data$sample == sample,], clinical_data = MSK_clinical_data[MSK_clinical_data$sample == sample,])
#' # Run INCOMMON classification
#' x = classify(x = x, priors = pcawg_priors, entropy_cutoff = NULL, rho = 0.01)
#' # An S3 method can be used to report to screen what is in the object
#' print(x)
#' @importFrom dplyr filter mutate rename select %>%

classify = function(
    x,
    k_max = 8,
    priors_m_k = pcawg_priors,
    priors_x = pcawg_priors_x,
    purity_error = 0.05,
    num_cores = NULL,
    iter_warmup = 500,
    iter_sampling = 1000,
    dump = FALSE,
    dump_file = NULL
)
  {

  stopifnot(inherits(x, "INCOMMON"))

  check_input(x)

  cli::cli_alert_info("Performing classification")

  classify_sample = function(x, sample){

    cli::cli_h1(
      "INCOMMON inference of copy number and mutation multiplicity for sample {.field {sample}}"
    )
    cat("\n")

    x = subset_sample(x = x, sample = sample)

    M = input(x) %>% nrow()
    N = input(x)$DP
    n = input(x)$NV
    purity_mean = purity(x = x, sample = sample)

    priors_m_k_sample = get_sample_priors(x = x, N_mutations = M, priors = priors_m_k, k_max = k_max)

    data = list(
      M = M,
      N = N,
      n = n,
      k_max = k_max,
      purity_mean = purity_mean,
      purity_error = purity_error,
      alpha_x = priors_x$alpha_x,
      beta_x = priors_x$beta_x,
      k_m_prior = priors_m_k_sample
    )

    model = get_stan_model()

    fit = model$sample(
      data = data,
      seed = 1992,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      parallel_chains = num_cores,
    )

    attach_fit_results(x = x, fit = fit)
  }

  lapply(samples(x), function(s){
    out = classify_sample(x = x, sample = s)

    if(!('classification' %in% names(x))) x$classification = list()
    if(!('fit' %in% names(x$classification))) x$classification$fit = dplyr::tibble(NULL)

    x$classification$fit <<- rbind(x$classification$fit, out)

    if(dump){
      x$classification$parameters = dplyr::tibble(
        k_max,
        purity_error,
        stan_iter_warmup = iter_warmup,
        stan_stan_iter_sampling = iter_sampling
      )

      x$classification$priors_m_k = priors_m_k
      x$classification$priors_x = priors_x

      if(is.null(dump_file)) dump_file = './dump.rds'

      saveRDS(object = x, file = dump_file)
    }

  })

  x$classification$parameters = dplyr::tibble(
    k_max,
    purity_error,
    stan_iter_warmup = iter_warmup,
    stan_stan_iter_sampling = iter_sampling
  )

  x$classification$priors_m_k = priors_m_k
  x$classification$priors_x = priors_x

    cli::cli_alert_info('There are: ')
    for (map_class in c('m=1', '1<m<k', 'm=k')) {
      N  = classification(x) %>% dplyr::filter(map_class == !!map_class) %>% nrow()
      cli::cli_bullets(c("*" = paste0("N = ", N, ' mutations (', map_class, ')')))
    }

    mean_ent = classification(x) %>% dplyr::pull(entropy) %>% mean()
    min_ent = classification(x) %>% dplyr::pull(entropy) %>% min()
    max_ent = classification(x) %>% dplyr::pull(entropy) %>% max()

    cli::cli_alert_info(
      'The mean classification entropy is {.field {round(mean_ent, 2)}} (min: {.field {round(min_ent, 2)}}, max: {.field {round(max_ent, 2)}})'
      )

  return(x)
}


