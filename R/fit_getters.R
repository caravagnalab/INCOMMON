get_fit_sample = function(x, sample){
  what = x$classification$fit %>% filter(sample == !!sample)
  out = what$fit[[1]]
  return(out)
}

sample = 'CPCT02010043T'

get_z_km = function(x, sample){
  what = get_fit_sample(x = x, sample = sample)
  x = subset_sample(x = x, sample_list = sample)
  z_km = what$summary(variables = 'z_km')  
  out_table = dplyr::tibble(NULL)
  for(i in 1:nrow(input(x))){
    idx = 1
    for(k in 1:k_max){
      for(m in 1:k){
        out_table = rbind(
          out_table,
          z_km %>% 
            filter(grepl(paste0('z_km\\[',i,',',idx,'\\]'), variable)) %>% 
            dplyr::mutate(k = k, m = m) %>% 
            dplyr::bind_cols(
              input(x)[i,] %>% dplyr::select(sample, chr, from, to, ref, alt))
        )
        idx = idx + 1
      }
    }
  }

  out_table %>% 
      dplyr::mutate(variable = 'z_km') %>% 
      dplyr::select(sample, chr, from, to, ref, alt, variable, k, m, dplyr::everything())
}

get_map_z_km = function(x, sample){
  get_z_km(x = x, sample = sample) %>% 
    dplyr::group_by(sample, chr, from, to, ref, alt) %>% 
    dplyr::arrange(desc(mean), .by_group = TRUE) %>% 
    dplyr::slice_head(n = 1)
}

get_map_purity = function(x, sample){
  what = get_fit_sample(x = x, sample = sample)
  purity = what$summary(variables = 'purity') %>% 
    dplyr::mutate(sample = sample)
  return(purity)
}


classification = function(x){
  samples = input(x) %>% dplyr::pull(sample) %>% unique()
  lapply(samples, function(s){
    get_map_z_km(x = x, sample = s) %>% 
      dplyr::rename(z_km = mean) %>% 
      dplyr::select(sample, chr, from, to, ref, alt, k, m, z_km) %>% 
      dplyr::full_join(
        input(x), 
        by = c('sample','chr', 'from', 'to', 'ref', 'alt')
        ) %>% 
      dplyr::full_join(
        get_map_purity(x, sample = s) %>% 
          dplyr::rename(purity_fit = mean) %>% 
          dplyr::select(sample, purity_fit)
      )
  }) %>% do.call(rbind, .)
}


classification(x) %>% View()
