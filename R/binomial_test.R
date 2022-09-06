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
                         karyotypes,
                         assign_extremes
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
        ploidy = Major+minor,
        multiplicity = p,
        karyotype = paste0(Major, ":", minor),
        label = paste0(Major+minor,'N (Mutated: ', p,"N)"),
        peak = peaks[p]
      )) %>%
      Reduce(f = bind_rows)
  }
  
  cut = function(x, cutoff)
  {
   
    ## Normalize likelihoods by maximum and compute corresponding cut-offs 
    x = x %>% 
      group_by(label) %>% 
      summarise(
        density,
        NV,
        ploidy,
        multiplicity,
        label,
        peak,
        cutoff = max(density) * cutoff
      ) %>% 
      unique() %>% 
      ungroup()
    
    ## Compute uncertainty (relative likelihood) for each class
    x = x %>% 
      group_by(NV) %>% 
      summarise(
        density,
        ploidy,
        multiplicity,
        peak,
        cutoff,
        label,
        uncertainty = 1 - max(density) / sum(density)
      ) %>% 
      ungroup()
    
    ## Add out-of-sample label for densities lower than cut-offs
    x = x %>% 
      mutate(label = ifelse(density < cutoff, "out of sample", label))
    x
  }
  
  
  dataset = lapply(karyotypes, function(k) {
    alleles = strsplit(k, split = ":")[[1]] %>% as.integer()
    db(alleles[1],alleles[2])
  }) %>% bind_rows() %>% cut(cutoff)
  
  tested = dataset %>%
    maximise() %>%
    dplyr::filter(NV == test)

    # Optionally assign NVs lower than minimum (larger than maximum) accepted to 
  # Subclonal/Higher Ploidy (Higher Ploidy and Multiplicity)
  if(tested$label == "out of sample" & assign_extremes){
    min_in_sample = dataset %>%
      dplyr::filter(label != "out of sample") %>%
      dplyr::filter(NV == min(NV)) %>%
      pull(NV)
    
    max_in_sample = dataset %>%
      dplyr::filter(label != "out of sample") %>%
      dplyr::filter(NV == max(NV)) %>%
      pull(NV)
    
    # tested$label = ifelse(tested$NV < min_in_sample | tested$NV > max_in_sample,
    #                       paste0(tested$ploidy,'N (Mutated: ', tested$multiplicity,"N)"),
    #                       tested$label)
    # tested$label = ifelse(tested$NV < min_in_sample | tested$NV > max_in_sample,
    #                       paste0(tested$ploidy,'N (Mutated: ', tested$multiplicity,"N)"),
    #                       tested$label)
    if(tested$NV < min_in_sample) tested$label = "SUB/HPLM"
    if(tested$NV > max_in_sample) tested$label = "HPHM"
  }
  
  
  return(tibble(ploidy = tested$ploidy,
                multiplicity = tested$multiplicity,
                uncertainty = tested$uncertainty,
                label = tested$label,
                density = list(dataset)))
}
  