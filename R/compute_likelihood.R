#' Compute model likelihood
#'
#' @param NV Number of reads with the variant.
#' @param DP Sequencing coverage of the mutated genome site.
#' @param prob Success probability (expected VAF).
#' @param rho The over-dispersion parameter.
#' @return A vector of probability densities (from NV = 1 to NV = DP).
#' @export
#' @examples
#'compute_likelihood(
#' NV = 170,
#' DP = 200,
#' prob = 0.5,
#' rho = 0.01
#')
compute_likelihood = function(NV, DP, prob, rho) {
  density = VGAM::dbetabinom(
    x = 1:DP,
    size = DP,
    prob = prob,
    rho = rho
  )
  return(density)
}
