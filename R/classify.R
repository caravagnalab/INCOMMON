#' Classify mutations using a Beta-Binomial model-based test.
#'
#' @param x An object of class \code{'INCOMMON'} generated with function `init`.
#' @param k_max The maximum value of total copy number to be included in the model.
#' @param priors_k_m A dplyr::tibble or data frame with columns `gene`, `gene_role`, `tumor_type`, `ploidy`, `multiplicity`
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
    priors_k_m = pcawg_priors,
    priors_x = pcawg_priors_x,
    purity_error = 0.05,
    num_cores = NULL,
    iter_warmup = 500,
    iter_sampling = 1000,
    num_chains = 4,
    dump = FALSE,
    dump_file = NULL,
    stan_fit_dump = TRUE,
    stan_fit_dir = './results/stan_fits'
)
  {

  stopifnot(inherits(x, "INCOMMON"))

  check_input(x)

  cli::cli_alert_info("Performing classification")

  classify_sample = function(x, sample, k_max){

    cli::cli_h1(
      "INCOMMON inference of copy number and mutation multiplicity for sample {.field {sample}}"
    )
    cat("\n")

    x = subset_sample(x = x, sample = sample)

    ## Prepare input

    M = input(x) %>% nrow()
    N = input(x)$DP
    n = input(x)$NV
    purity_mean = purity(x = x, sample = sample)
    tumor_type = input(x) %>% dplyr::pull(tumor_type) %>% unique()

    priors_k_m_sample = get_stan_input_priors(x = x, N_mutations = M, priors = priors_k_m, k_max = k_max)
    priors_x_sample = priors_x %>% filter(tumor_type == !!tumor_type)
    if(nrow(priors_x_sample)==0){
      priors_x_sample = priors_x %>% filter(tumor_type == 'PANCA')
    }

    data = list(
      M = M,
      N = N,
      n = n,
      k_max = k_max,
      purity_mean = purity_mean,
      purity_error = purity_error,
      alpha_x = priors_x_sample$alpha_x,
      beta_x = priors_x_sample$beta_x,
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

    ## Draw samples from priors and posteriors

    k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
      dplyr::as_tibble() %>%
      dplyr::filter(m <= k) %>% arrange(k, m)

    N_rep = fit$draws(variables = 'N_rep') %>% array(dim = c(num_chains * iter_sampling, M))
    n_rep = fit$draws(variables = 'n_rep') %>% array(dim = c(num_chains * iter_sampling, M))

    x_rep = fit$draws(variables = 'x') %>% array()
    x_prior_rep = fit$draws(variables = 'x_rep') %>% array()

    purity_rep = fit$draws(variables = 'purity') %>% array()
    purity_prior_rep = fit$draws(variables = 'purity_rep') %>% array()

    km_idx = fit$draws(variables = 'km_idx') %>% array(dim = c(num_chains * iter_sampling, M))
    km_rep = lapply(1:M, function(i){
      lapply(1:length(km_idx[,i]), function(j){
        k_m_table[km_idx[j,i],]
      }) %>% do.call(rbind, .)
    })

    z_km = fit$summary(variables = 'z_km')
    z_km = lapply(1:M, function(i){
      z_km %>%
        dplyr::filter(grepl(paste0('z_km\\[',i), variable)) %>%
        dplyr::bind_cols(k_m_table) %>%
        dplyr::select(k, m, median) %>%
        dplyr::rename(z_km = median)
    })

    ## MAP estimates

    km_map = lapply(km_rep, function(x){
      what = x %>% table() %>% dplyr::as_tibble() %>% dplyr::mutate(k = as.integer(k), m = as.integer(m))
      dplyr::arrange(what, dplyr::desc(n)) %>% dplyr::slice_head(n=1) %>% dplyr::select(-n)
    })

    x_map = median(x_rep)
    purity_map = median(purity_rep)

    ## Bayesian p-value tests

    bayes_p_x = bayesian_p_value(posterior_rep = x_rep, prior_rep = x_prior_rep, prior_type = 'gamma')
    bayes_p_purity = bayesian_p_value(posterior_rep = purity_rep, prior_rep = purity_prior_rep, prior_type = 'beta')

    ## Generate report

    ppois = plot_poisson_model(x =  x, sample = sample, N_rep = N_rep, km_rep = km_rep, km_map = km_map, purity_map = purity_map, x_map = x_map)

    px = plot_x_check(posterior_rep = x_rep, prior_rep = x_prior_rep, bayes_p = bayes_p_x)
    ppi = plot_purity_check(posterior_rep = purity_rep, prior_rep = purity_prior_rep, bayes_p = bayes_p_purity)
    p_x_pi = patchwork::wrap_plots(px, ppi, design = 'A\nB', guides = 'collect')&ggplot2::theme(legend.position = 'bottom')

    pbin = plot_binomial_model(x = x, n_rep = n_rep, km_rep = km_rep)

    pkmprior = plot_prior_k_m(priors_k_m = priors_k_m, x = x, k_max = k_max)
    pkmpost = plot_posterior_k_m(x = x, M = , z_km = z_km)

    plot_report = patchwork::wrap_plots(
      ppois, p_x_pi, pkmprior, pkmpost, pbin,
      design = 'AAABBB\nCCDDEE\nCCDDEE'
      ) +
      patchwork::plot_annotation(
        title = paste0(sample),
        subtitle = paste0(tumor_type, '; M = ', M, ' mutations'))


    if(stan_fit_dump){

      if (!dir.exists(stan_fit_dir)) {
        dir.create(stan_fit_dir, recursive = TRUE)
      }

      fit$save_object(file = paste0(stan_fit_dir, '/', sample, '.rds'))
    }

    km_map %>% do.call(rbind, .) %>% dplyr::mutate

    tibble(
      purity_map,
      bayes_p_purity,
      x_map,
      bayes_p_x,
      z_km
    ) %>%
      dplyr::bind_cols(
        km_map %>% do.call(rbind, .) %>% dplyr::rename(map_k = k, map_m = m)
      )
}

  lapply(samples(x), function(s){
    fit = classify_sample(x = x, sample = s, k_max = k_max)

    # idx = 1
    # class_table = tibble(NULL)
    # for(k in 1:k_max){
    #   for(m in 1:k){
    #     class_table = class_table %>%
    #       rbind(tibble(idx = idx, k = k,m = m))
    #     idx = idx+1
    #   }
    # }

    # out %>%
    #   dplyr::rowwise() %>%
    #   dplyr::mutate(
    #     idx = strsplit(variable, split = ',')[[1]][2] %>%
    #       gsub('\\]', '', .) %>% as.integer()
    #   ) %>% full_join(class_table, by = 'idx') %>%
    #   dplyr::select(sample, k, m, dplyr::everything())


    if(!('classification' %in% names(x))) x$classification = list()
    if(!('fit' %in% names(x$classification))) x$classification$fit = list()

    # x$classification$fit <<- rbind(x$classification$fit, dplyr::tibble(sample = s, fit = list(fit)))
    x$classification$fit[[s]] <<- fit

    x$classification$parameters <<- dplyr::tibble(
      k_max,
      purity_error,
      stan_iter_warmup = iter_warmup,
      stan_iter_sampling = iter_sampling,
      num_chains,
      dump,
      dump_file,
      stan_fit_dump,
      stan_fit_dir
    )

    x$classification$priors_k_m <<- priors_k_m
    x$classification$priors_x <<- priors_x

    if(dump){

      if(is.null(dump_file)) dump_file = './dump.rds'

      saveRDS(object = x, file = dump_file)
    }

  })

    # cli::cli_alert_info('There are: ')
    # for (map_class in c('m=1', '1<m<k', 'm=k')) {
    #   N  = classification(x) %>% dplyr::filter(map_class == !!map_class) %>% nrow()
    #   cli::cli_bullets(c("*" = paste0("N = ", N, ' mutations (', map_class, ')')))
    # }

    # mean_ent = classification(x) %>% dplyr::pull(entropy) %>% mean()
    # min_ent = classification(x) %>% dplyr::pull(entropy) %>% min()
    # max_ent = classification(x) %>% dplyr::pull(entropy) %>% max()

    # cli::cli_alert_info(
    #   'The mean classification entropy is {.field {round(mean_ent, 2)}} (min: {.field {round(min_ent, 2)}}, max: {.field {round(max_ent, 2)}})'
    #   )

  return(x)
}


