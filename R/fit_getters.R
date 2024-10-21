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
  samples = input(x) %>% dplyr::pull(sample) %>% unique()
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
