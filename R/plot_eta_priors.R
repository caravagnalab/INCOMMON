#' Visualise the prior distribution over the rate of reads per chromosome copy.
#'
#' @param priors_eta An object of class INCOMMON.
#' @export
#' @examples
#' # Load the default prior on eta, estimated from the MKS-MET data.
#' data(priors_eta)
#' # Plot classification results for a specific sample
#' plot_eta_prior(priors_eta = priors_eta)
#' @importFrom dplyr filter mutate group_by reframe select %>%
#' @importFrom stats var
plot_eta_prior = function(priors_eta){
  toplot = lapply(1:nrow(priors_eta), function(i){
    tibble(
      tumor_type = priors_eta[i,]$tumor_type,
      cohort = priors_eta[i,]$cohort,
      x = rgamma(n = 1:10000, shape = priors_eta[i,]$alpha_eta, rate = priors_eta[i,]$beta_eta)
    )
  }) %>% do.call(rbind, .)

  toplot_dp = x$input %>% select(DP, tumor_type)

  rbind(
    toplot %>% mutate(what = 'Read count rate per copy (Prior)'),
    toplot_dp %>% mutate(what = 'Sequencing Depth (Oberved)') %>% rename(x = DP)
  ) %>%
    filter(
      tumor_type %in% priors_eta$tumor_type
    ) %>%
    mutate(
      tumor_type = factor(tumor_type)
    ) %>%
    mutate(tumor_type = relevel(tumor_type, ref = 'PANCA')) %>%
    ggplot2::ggplot(ggplot2::aes(x = x, fill = what))+
    ggridges::geom_density_ridges(ggplot2::aes(y = tumor_type), alpha = .5, scale = .9)+
    ggplot2::scale_fill_manual(values = c('indianred', 'steelblue'))+
    my_ggplot_theme()+
    ggplot2::labs(
      x = '',
      y = 'density',
      fill = ''
    )+
    ggplot2::xlim(0,1500)

}
