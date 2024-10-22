get_z_km = function(x, sample){
  what = x$classification$fit[[sample]]
  z_km = what$z_km %>% colMeans()
  k_max = x$classification$parameters$k_max

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

get_map_z_km = function(x, sample){
  get_z_km(x = x, sample = sample) %>%
    dplyr::group_by(i) %>%
    dplyr::arrange(desc(z_km), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-i)
}

get_map_purity = function(x, sample){
  x$classification$fit[[sample]]$purity_fit %>% median()
}

get_map_x = function(x, sample){
  x$classification$fit[[sample]]$x %>% median()
}


classification = function(x){
  samples = names(x$classification$fit)
  x = subset_sample(x = x, sample_list = samples)
  lapply(samples, function(s){
    what = get_map_z_km(x = x, sample = s)
    # what$sample = s
    what$purity_fit = get_map_purity(x = x, sample = s)
    what$x_fit = get_map_x(x = x, sample = s)
    what
  }) %>% do.call(rbind, .) %>%
    dplyr::bind_cols(x$input) %>%
    dplyr::select(sample, purity, chr, from, to , ref, alt, gene, gene_role, NV, DP, VAF, k, m, z_km, purity_fit, x_fit, dplyr::everything())
}

get_N_rep = function(x, sample){
  x$classification$fit[[sample]]$N_rep
}

get_n_rep = function(x, sample){
  x$classification$fit[[sample]]$n_rep
}


get_N_rep_ci = function(x, sample) {
  N_rep = get_N_rep(x = x, sample = sample) %>% as.vector()
  dplyr::tibble(
    N = N_rep,
    k = get_k_m_rep(x = x, sample = sample)$k
  ) %>%
    dplyr::group_by(k) %>%
    dplyr::reframe(mean = mean(N), q5 = quantile(N, probs = .05), q95 = quantile(N, probs = .95))
}

get_k_m_rep = function(x, sample){

  k_max = x$classification$parameters$k_max

  z_km = x$classification$fit[[sample]]$z_km

  k_m_table = expand.grid(k = 1:k_max, m = 1:k_max) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(m <= k) %>% arrange(k, m)

  km_idx = apply(z_km, c(1, 2), which.max)

  dplyr::tibble(
    k = k_m_table$k[as.vector(km_idx)],
    m = k_m_table$m[as.vector(km_idx)]
  )

}



