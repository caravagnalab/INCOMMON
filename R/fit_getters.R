get_map_z_km = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)
  what = x$classification$fit[[sample]]
  M = nrow(input(x))

  k_max = x$classification$parameters$k_max

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% arrange(k, m)

  lapply(1:M, function(i){
    z_km = what %>%
      dplyr::filter(grepl(paste0('z_km\\[',i), variable)) %>%
      dplyr::arrange(dplyr::desc(median)) %>%
      dplyr::slice_head(n = 1)

    map_km = what %>%
      dplyr::filter(grepl(paste0('km_idx\\[',i,'\\]'), variable)) %>%
      dplyr::pull(median) %>% as.integer()

    tibble(
      sample = sample,
      k = k_m_table[map_km,]$k,
      m = k_m_table[map_km,]$m,
      z_km = z_km$median
    )
  }) %>% do.call(rbind, .)
}

get_map_purity = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)
  what = x$classification$fit[[sample]]

  purity_fit = what %>% dplyr::filter(grepl('purity', variable)) %>% dplyr::pull(median)

  tibble(
    sample = sample,
    purity_fit
  )
}

get_map_count_rate = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)
  what = x$classification$fit[[sample]]

  x_fit = what %>% dplyr::filter(grepl('^x$', variable)) %>% dplyr::pull(median)

  tibble(
    sample = sample,
    x_fit
  )
}


classification = function(x){
  samples = names(x$classification$fit)
  x = subset_sample(x = x, sample_list = samples)

  lapply(samples, function(s){
    list(
      get_map_z_km(x = x, sample = s),
      get_map_purity(x = x, sample = s),
      get_map_count_rate(x = x, sample = s)
    ) %>% Reduce(function(x,y) full_join(x,y, by = 'sample'), .)
      }) %>% do.call(rbind, .) %>%
    dplyr::select(-sample) %>%
    dplyr::bind_cols(x$input) %>%
    dplyr::select(sample, purity, chr, from, to , ref, alt, gene, gene_role, NV, DP, VAF, k, m, z_km, purity_fit, x_fit, dplyr::everything())
}


get_z_km_draws = function(x, sample){

  x = subset_sample(x = x, sample_list = sample)

  stopifnot(x$classification$parameters$stan_fit_dump)

  dir = x$classification$parameters$stan_fit_dir
  num_chains = x$classification$parameters$num_chains
  iter_sampling = x$classification$parameters$stan_iter_sampling
  k_max = x$classification$parameters$k_max

  M = classification(x) %>% nrow()

  fit = readRDS(paste0(dir, '/', sample, '.rds'))

  km_idx = fit$draws(variables = 'km_idx')
  km_idx = array(km_idx, dim = c(num_chains * iter_sampling, M))

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% arrange(k, m)

  z_km = lapply(1:nrow(z_km), function(i){
    k_m_table$z_km = z_km[i,]
    k_m_table$i = i
    k_m_table
  }) %>% do.call(rbind, .)

 z_km
}
#

get_N_draws = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)

  stopifnot(x$classification$parameters$stan_fit_dump)

  dir = x$classification$parameters$stan_fit_dir
  num_chains = x$classification$parameters$num_chains
  iter_sampling = x$classification$parameters$stan_iter_sampling

  M = classification(x) %>% nrow()

  fit = readRDS(paste0(dir, '/', sample, '.rds'))

  N_rep = fit$draws(variables = 'N_rep')
  N_rep = array(N_rep, dim = c(num_chains * iter_sampling, M))
}

get_n_draws = function(x, sample){
  x = subset_sample(x = x, sample_list = sample)

  stopifnot(x$classification$parameters$stan_fit_dump)

  dir = x$classification$parameters$stan_fit_dir
  num_chains = x$classification$parameters$num_chains
  iter_sampling = x$classification$parameters$stan_iter_sampling

  M = classification(x) %>% nrow()

  fit = readRDS(paste0(dir, '/', sample, '.rds'))

  n_rep = fit$draws(variables = 'n_rep')
  n_rep = array(N_rep, dim = c(num_chains * iter_sampling, M))
}


get_N_rep_ci = function(x, sample) {
  N_rep = get_N_draws(x = x, sample = sample) %>% as.vector()
  dplyr::tibble(
    N = N_rep,
    k = get_k_m_draws(x = x, sample = sample)$k
  ) %>%
    dplyr::group_by(k) %>%
    dplyr::reframe(mean = mean(N), q5 = quantile(N, probs = .05), q95 = quantile(N, probs = .95))
}
#
get_k_m_draws = function(x, sample){

  k_max = x$classification$parameters$k_max

  dir = x$classification$parameters$stan_fit_dir
  # num_chains = x$classification$parameters$num_chains
  # iter_sampling = x$classification$parameters$stan_iter_sampling
  M = classification(x) %>% nrow()

  fit = readRDS(paste0(dir, '/', sample, '.rds'))

  km_idx = fit$draws(variables = 'km_idx')
  km_idx = array(km_idx, dim = c(num_chains * iter_sampling, M))

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% arrange(k, m)

  y = lapply(1:M, function(i){
    lapply(1:length(km_idx[,i]), function(j){
        k_m_table[km_idx[j,i],]
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

  y

}
