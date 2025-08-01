#' @param x An object of class INCOMMON.
#' @param priors_k_m Pre-computed priors for gene muation total copy number and multipliity.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Plot classification results for a specific sample
#' plot_prior(x = MSK_classified, purity_error = 0.05)
#' @importFrom dplyr filter mutate group_by reframe select %>%
#' @importFrom stats var
#'
compute_eta_prior = function(x, priors_k_m){

  priors_k_m = priors_k_m %>%
    dplyr::group_by(gene, tumor_type) %>%
    dplyr::mutate(p = n / sum(n))

  x_data = x$input %>%
    dplyr::select(sample, tumor_type, gene, purity, DP) %>%
    unique() %>%
    dplyr::left_join(priors_k_m, by = c('gene', 'tumor_type')) %>%
    dplyr::mutate(x = DP*p/(2*(1-purity)+k*purity)) %>%
    dplyr::group_by(tumor_type, gene, sample) %>%
    dplyr::reframe(x = sum(x, na.rm = T)) %>%
    dplyr::filter(!is.na(x)) %>% filter(x > 0) %>%
    unique() %>%
    dplyr::group_by(sample) %>%
    dplyr::reframe(x = mean(x), tumor_type) %>%
    unique()

  prior_x_ts = x_data %>%
    dplyr::group_by(tumor_type) %>%
    dplyr::reframe(mean_eta = mean(x), var_eta = stats::var(x), N = length(unique(sample))) %>%
    dplyr::filter(N >= 50) %>%
    dplyr::mutate(
      alpha_eta = (mean_eta**2)/var_eta,
      beta_eta = mean_eta/var_eta
    ) %>%
    dplyr::filter(!is.infinite(alpha_eta))

  prior_x_pan = x_data %>%
    dplyr::mutate(tumor_type = 'PANCA') %>%
    unique() %>%
    dplyr::group_by(tumor_type) %>%
    dplyr::reframe(mean_eta = mean(x), var_eta = stats::var(x), N = length(unique(sample))) %>%
    dplyr::mutate(
      alpha_eta = (mean_eta**2)/var_eta,
      beta_eta = mean_eta/var_eta
    ) %>%
    dplyr::filter(!is.infinite(alpha_eta))

  prior_eta = rbind(prior_x_ts, prior_x_pan)

  return(prior_eta)
}
