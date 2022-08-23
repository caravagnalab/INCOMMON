#' Compute likelihoods and uncertainty..
#'
#' @param test Number of reads with the variant.
#' @param DP Sequencing coverage of the mutated genome site.
#' @param purity Purity of the sample.
#' @param cutoff Likelihood cut-off for class assignment.
#' @param model Model used for the classification task, either "Binomial" (no over-dispersion) 
#' or "Beta-Binomial" (over-dispersion included), that will be used as the expected 
#' distribution for the number of reads with variant at fixed coverage and purity.
#' @param rho If "Beta-Binomial" model is selected, this parameter tunes the over-dispersion
#' of the expected distribution used for the classification.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @return A tibble including ploidy, multiplicity, proportion of wt alleles, 
#' uncertainty  of the class assignment, model density.
#' @export
#' @importFrom dplyr filter mutate rename select %>% 
#' @examples
#' binomial_test(test = 220,
#' DP = 500,
#' purity = 0.98,
#' cutoff = 0.9,
#' model = "Beta-Binomial",
#' rho = 0.01,
#' karyotypes = c("1:0","1:1","2:0","2:1","2:2")
#')

binomial_test = function(test,
                         DP,
                         purity,
                         cutoff,
                         model,
                         rho = 0.01,
                         karyotypes
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
        label = paste0(Major, ":", minor, ' ', p),
        peak = peaks[p]
      )) %>%
      Reduce(f = bind_rows)
  }
  
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
  
  dataset = lapply(karyotypes, function(k) {
    alleles = strsplit(k, split = ":")[[1]] %>% as.integer()
    db(alleles[1],alleles[2])
  }) %>% bind_rows() %>% cut(cutoff)
  
  normll = dataset %>% 
    dplyr::mutate(state = paste0(Major, ":", minor, " ", multiplicity)) %>% 
    group_by(state) %>%
    dplyr::summarise(NV,density = density/max(density)) %>%
    dplyr::filter(NV==test) %>% 
    dplyr::arrange(desc(density))
    
  
  llratio = normll[2,]$density/normll[1,]$density
  
  class_of = dataset %>%
    maximise() %>%
    dplyr::filter(NV == test) %>% 
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
                uncertainty = llratio,
                density = list(dataset)))
}
  