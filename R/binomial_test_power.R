# test - scalare della mutazione da testare
# DP - coverage della mutazione da testare
# purity - del campione
# cutoff - likelihood minima da considerare per il test
# normalise - se ogni likelihood viene ri-scalato sul suo massimo (nota, i valori
# sono dapo tra 0/1 ma NON sono probability, sono re-scaled likelihoods)
binomial_test_power = function(test, 
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
  
  # Maxima
  maximise = function(x)
  {
    x %>%
      group_by(NV) %>%
      filter(density == max(density)) %>%
      ungroup()
  }
  
  # Re-scaled likelihoods (makes them comparable)
  scaling = function(x, normalise)
  {
    if (!normalise)
      return(x)
    
    x %>%
      group_by(karyotype, multiplicity) %>%
      mutate(density = density / max(density)) %>%
      ungroup()
  }
  
  # Cutoff
  cut = function(x, cutoff)
  {
    x %>%
      mutate(label = ifelse(density < cutoff, 'out of sample', label))
  }
  
  
  colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
  names_colors = expand.grid(names(colors), 1:2) %>% apply(1, paste, collapse = ' ')
  colors = sapply(names_colors, function(n)
    colors[strsplit(n, ' ')[[1]][1]])
  names(colors) = names_colors
  
  colors = c(colors, `out_of_sample` = 'gray')
  
  dataset = bind_rows(db(1, 0),
                      db(1, 1),
                      db(2, 0),
                      db(2, 1),
                      db(2, 2)) %>% scaling(normalise = TRUE) %>% cut(cutoff)
  
  class_of =  dataset %>%
    maximise() %>%
    filter(NV == test) %>%
    pull(label) %>% 
    unique() %>% 
    paste(collapse = ', ')
  
  # plot_test =
  plot_test =
    dataset %>%
    ggplot() +
    geom_line(aes(x = NV, y = density), color = 'gainsboro',  size = .3) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    geom_point(data = dataset %>% maximise(),
               aes(x = NV,
                   y = density,
                   color = label),
               size = 1) +
    geom_line(data = dataset %>% maximise(),
              aes(x = NV,
                  y = density,
                  color = label),
              size = .3) +
    labs(title = paste("Class:", class_of),
         subtitle = paste("Purity", purity, ' - cutoff', cutoff)) +
    guides(color = 'none') +
    geom_hline(
      yintercept = cutoff,
      linetype = 'dashed',
      color = 'gray',
      size = .5
    ) +
    geom_vline(
      xintercept = test,
      color = ifelse('out_of_sample' %in% class_of, 'indianred3', 'forestgreen'),
      linetype = 'dashed',
      size = .5
    )
    # facet_wrap(~karyotype)
    # facet_grid(cols = vars(multiplicity), rows = vars(karyotype))
  tibble(
    class = class_of,
    likelihood = dataset %>%
      maximise() %>%
      filter(NV == test) %>%
      pull(density) %>% unique(),
    plot = list(plot_test)
  )
}

# binomial_test_power(test = 320, DP =  500, cutoff = .1,
#                     purity =  1)
# binomial_test_power(80, 240, .97, .04, normalise = FALSE)
# binomial_test_power(80, 240, .7, .9)