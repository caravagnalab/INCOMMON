#' Compute model posterior and entropy
#'
#' @param NV Number of reads with the variant.
#' @param DP Sequencing coverage of the mutated genome site.
#' @param gene Gene name or symbol.
#' @param priors Prior distribution.
#' @param purity Purity of the sample.
#' @param tumor_type Tumor type of the sample.
#' @param entropy_cutoff Cut-off on entropy for Tier-1/Tier-2 distinction.
#' @param rho The over-dispersion parameter.
#' @param karyotypes Karyotypes to be included among the possible classes.
#' @return A table including ploidy, multiplicity, posterior probability,
#' and classification entropy.
#' @export
#' @importFrom dplyr filter mutate rename select %>%
#' @examples
#'compute_posterior(
#' NV = 170,
#' DP = 200,
#' gene = 'TP53',
#' priors = pcawg_priors,
#' tumor_type = 'PAAD',
#' purity = 0.9,
#' entropy_cutoff = 0.2,
#' rho = 0.01,
#' karyotypes = c("1:0", "1:1", "2:0", "2:1", "2:2")
#')

compute_posterior = function(NV,
                         DP,
                         gene,
                         priors = NULL,
                         tumor_type,
                         purity,
                         entropy_cutoff,
                         rho = 0.01,
                         karyotypes,
                         silent = FALSE)
{
  NV_x = 1:DP

  # Density
  db = function(Major, minor, prior, gene)
  {

    expected_peaks = CNAqc:::expected_vaf_peak(Major, minor, purity)$peak

    lapply(expected_peaks %>% seq_along(), function(p) {

      # Classification label
      label = paste0(Major + minor, 'N (Mutated: ', p, "N)")

      if(is.data.frame(prior)){
        if(!(label %in% prior$label) & !silent) cli::cli_alert_danger("Incomplete prior distribution!")
        stopifnot(label %in% prior$label)
        prior = prior %>% dplyr::filter(label == !!label) %>% dplyr::pull(p)
      }

      # Expected VAF peak of mixture component
      expected_peak = expected_peaks[p]

      # Beta-Binomial likelihood distribution
      likelihood = compute_likelihood(NV = NV_x, DP = DP, prob = expected_peak, rho = rho)

      # # Posterior distribution
      # if (is.null(priors)){
      #   prior = 1
      # } else {
      #   prior = get_prior(priors, gene, tumor_type, label)
      # }

      posterior = prior*likelihood

      # Output
      out = data.frame(
        value = posterior,
        NV = NV_x,
        Major = Major,
        minor = minor,
        ploidy = Major + minor,
        multiplicity = p,
        karyotype = paste0(Major, ":", minor),
        label = label,
        peak = expected_peak
      )

      out

    }) %>%
      do.call(dplyr::bind_rows, .)
  }

  if(is.na(purity)){
    if(!silent) cli_alert_warning(text =
                        "With purity {.field {purity}} classification is not possible."
    )
    return(dplyr::tibble(ploidy = NA,
                  multiplicity = NA,
                  entropy = NA,
                  label = NA,
                  density = list(NULL)))
  }

  # Compute posterior distribution of 0 < NV <= DP for each component of the mixture (karyotype)
  # Prior distribution
  if (is.null(priors)){
    prior = 1
  } else {
    prior = get_prior(priors, gene, tumor_type, silent = silent)
  }
  posterior = lapply(karyotypes, function(k) {
    alleles = strsplit(k, split = ":")[[1]] %>% as.integer()
    db(Major = alleles[1],
       minor = alleles[2],
       prior =  prior,
       gene = gene)
  }) %>%
    dplyr::bind_rows()

  # Compute entropy
  posterior = posterior %>%
    dplyr::group_by(NV) %>%
    dplyr::reframe(value, NV, Major, minor, ploidy, multiplicity, karyotype, label, peak,
                   entropy = -sum(value*log2(max(value, .Machine$double.xmin)))) # entropy formula replacing zeros with minimum
                                                                                 # machine floating number

  # Apply entropy cutoff
  posterior = posterior %>%
    dplyr::mutate(state = dplyr::case_when(
      label %in% c('4N (Mutated: 1N)', '3N (Mutated: 1N)') | entropy > entropy_cutoff ~ 'Tier-2',
      entropy <= entropy_cutoff & label %in% c("2N (Mutated: 1N)") ~ "HMD",
      entropy <= entropy_cutoff & label %in% c("1N (Mutated: 1N)") ~ "LOH",
      entropy <= entropy_cutoff & label %in% c("2N (Mutated: 2N)") ~ "CNLOH",
      entropy <= entropy_cutoff & label %in% c("3N (Mutated: 2N)", "4N (Mutated: 2N)") ~ "AM"
    ))

  return(posterior)
}

