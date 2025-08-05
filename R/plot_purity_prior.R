#' Visualise the prior distribution on sample purity
#' @param x An object of class INCOMMON.
#' @param purity_error The variance of the Beta prior distribution on purity.
#' @param sample An identifier of the tumour sample.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_PAAD_output)
#' # Plot classification results for a specific sample
#' plot_purity_prior(x = MSK_PAAD_output, sample = "P-0000142", purity_error = 0.05)
#' @importFrom dplyr filter mutate rename select %>% tibble intersect
#' @importFrom stats rgamma
plot_purity_prior = function(x, sample, purity_error = 0.05){

  stopifnot(length(dplyr::intersect(samples(x), sample))>0)
  purity_mean = purity(x = x, sample = sample)
  alpha_pi = purity_mean * ((purity_mean * (1 - purity_mean) / purity_error) - 1)
  beta_pi = (1 - purity_mean) * (purity_mean * (1 - purity_mean) / purity_error)

  data = tibble(x = stats::rgamma(n = 100000, shape = alpha_pi, rate = beta_pi))

  toplot = dplyr::tibble(
    x = seq(0, 1, length.out = 1000),
    density = stats::dbeta(
      seq(0, 1, length.out = 1000),
      shape1 = alpha_pi,
      shape2 = beta_pi)
  )

  # Plot using ggplot2
  toplot %>%
    ggplot2::ggplot(ggplot2::aes(x = x, y = density)) +
    ggplot2::geom_line(color = 'steelblue', size = 1.2) +
    ggplot2::geom_vline(xintercept = purity_mean, linetype = 'longdash')+
    ggplot2::labs(
      title = paste0("Prior Purity Distribution"),
      subtitle = paste0("Beta Distribution (α = ", alpha_pi, ", β = ", beta_pi, "); Mean = ", purity_mean),
      x = bquote(pi),
      y = "Density"
    )+
    my_ggplot_theme()

}
