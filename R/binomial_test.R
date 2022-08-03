# test - scalare della mutazione da testare
# DP - coverage della mutazione da testare
# purity - del campione
# cutoff - likelihood minima da considerare per il test
# normalise - se ogni likelihood viene ri-scalato sul suo massimo (nota, i valori
# sono dapo tra 0/1 ma NON sono probability, sono re-scaled likelihoods)
binomial_test = function(test, 
                               DP, 
                               purity, 
                               cutoff, 
                               normalise = TRUE, 
                               model,
                               rho = 0.01
)
{
  NV_x = 1:DP
  
  # Density
  db = function(Major, minor)
  {
    peaks = CNAqc:::expected_vaf_peak(Major, minor, purity)$peak
    
    lapply(peaks %>% seq_along, function(p)
      data.frame(
        density = compute_density(NV_x, DP, prob = peaks[p], model, rho),
        NV = NV_x,
        Major = Major,
        minor = minor,
        multiplicity = p,
        karyotype = paste0(Major, ":", minor),
        label = paste0(Major, ":", minor, ' ', p)
      )) %>%
      Reduce(f = bind_rows)
  }
  
  # Re-scaled likelihoods (makes them comparable)
  scaling = function(x, normalise)
  {
    x %>%
      group_by(karyotype, multiplicity) %>%
      mutate(density = density / max(density)) %>%
      ungroup()
  }
  
  # Cutoff
  cut = function(x, cutoff)
  {
    cut_offs = x %>%  
                  group_by(label) %>% 
                  filter(density == max(density)) %>% 
                  ungroup() %>% 
                  mutate(cutoff = density*cutoff) %>% 
                  select(label, cutoff) 
    x %>% 
      left_join(cut_offs, by = "label") %>% 
      mutate(label = ifelse(density < cutoff, "out of sample", label))
  }
  
  dataset = bind_rows(db(1, 0),
                      db(1, 1),
                      db(2, 0),
                      db(2, 1),
                      db(2, 2)) %>% cut(cutoff)
  
  class_of =  dataset %>%
    maximise() %>%
    filter(NV == test) %>%
    pull(label) %>% 
    unique() %>% 
    paste(collapse = ', ')
  
  ploidy = NA
  multiplicity = NA
  if(class_of != "out of sample"){
      info = strsplit(class_of, ",")[[1]][1] %>% strsplit(" ")
      ploidy = strsplit(info[[1]][1],"\\:")[[1]] %>% as.integer() %>% sum()
      multiplicity = info[[1]][2] %>% as.integer()
  }
  
  return(tibble(ploidy = ploidy,
                multiplicity = multiplicity,
                wt = ploidy-multiplicity,
                density = list(dataset)))
}
  