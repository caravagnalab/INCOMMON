#' Compute model likelihood
#'
#' @param NV Number of reads with the variant.
#' @param DP Sequencing coverage of the mutated genome site.
#' @param prob Success probability (expected VAF).
#' @param rho The over-dispersion parameter.
#' @param purity_error Standard deviation of the Beta distribution over the
#' sample purity measure.
#' @param purity The estimated purity of the sample.
#' @return A vector of probability densities (from NV = 1 to NV = DP).
#' @export
#' @importFrom purrr reduce
#' @examples
#'compute_likelihood(
#' NV = 170,
#' DP = 200,
#' prob = 0.5,
#' rho = 0.01
#')
compute_likelihood = function(NV, DP, m, pl, prob, rho, purity, purity_error) {

  # Functions for the compund beta model
  expected_vaf = function(pu, m, pl){
    f = (m * pu)/(2*(1 - pu) + (pl * pu))
    return(f)
  }

  alpha_beta_pi <- function(purity_mean, purity_variance) {
    alpha = purity_mean * (((purity_mean * (1 - purity_mean)) / (purity_variance)) - 1)
    beta = (1 - purity_mean) * (((purity_mean * (1 - purity_mean)) / (purity_variance)) - 1)
    return(c(alpha, beta))
  }

  # Function to compute the alpha_f(p) and beta_f(p) for f | p, rho
  alpha_beta_f <- function(pu, rho, m, pl) {
    mean_f <- expected_vaf(pu = pu, m = m, pl = pl)
    precision <- (1 - rho) / rho
    alpha <- precision * mean_f
    beta <- precision * (1 - mean_f)
    return(c(alpha, beta))
  }

  # Function to evaluate the integrand for P(x) over f and p
  integrand_beta <- function(x, N, f, pu, m, pl, rho, purity, purity_error) {

    ab_f = alpha_beta_f(pu = pu, rho = rho, m = m, pl = pl) # fix over-dispersion coefficient and dependency of f on pi
    alpha_f <- ab_f[1]
    beta_f <- ab_f[2]

    ab_pi = alpha_beta_pi(purity_mean = purity, purity_variance = purity_error**2) # fix mean purity to input one and variance to fixed value
    alpha_pi <- ab_pi[1]
    beta_pi <- ab_pi[2]

    # Probability mass function of binomial and density of Beta for f
    out = dbinom(x, size = N, prob = f) * dbeta(f, alpha_f, beta_f) * dbeta(pu, alpha_pi, beta_pi)
    out[is.na(out)] = 0

    return(out)
  }

  numerical_integral = function(stride = 0.01,
                                NV,
                                DP,
                                m,
                                pl,
                                rho,
                                purity,
                                purity_error) {

    vals = expand.grid(purity = seq(0, 1 - 1e-16, stride),
                       freq = seq(0, 1, stride)) %>% dplyr::as_tibble()

    lapply(1:nrow(vals), function(i) {

      pu = vals[i, ]$purity
      f = vals[i, ]$freq

      integrand_beta(
        x = NV,
        N = DP,
        f = f,
        pu = pu,
        m = m,
        pl = pl,
        rho = rho,
        purity = purity,
        purity_error = purity_error
      ) * stride * stride
    }) %>% purrr::reduce(., `+`)
  }

  if(purity_error == 0){
    density = VGAM::dbetabinom(
      x = NV,
      size = DP,
      prob = prob,
      rho = rho
    )
  } else {

    max_purity_error = sqrt((purity)*(1 - purity))

    if(purity_error >= max_purity_error){
      cli_alert_warning(
        text ="Maximum purity error is {.field {round(max_purity_error, 2)}}, input {.field {purity_error}}")
      }

    density = numerical_integral(
      stride = 0.01,
      NV = NV,
      DP = DP,
      m = m,
      pl = pl,
      rho = rho,
      purity = purity,
      purity_error = purity_error
    )
  }

  return(density)
}
