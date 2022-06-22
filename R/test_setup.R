test_setup = function(coverage = 500,
                          purity = 1.0,
                          rho = 0.01,
                          alpha_level = 0.01,
                          model = 'binomial')
{
  # Range of NV values
  nvs = 1:coverage

  # Compute Binomial or Beta-Binomial probability for NV values in range
  log_p = NULL
  p = purity/2
  if (model == 'binomial')
  {
    log_p = sapply(nvs, dbinom, size = coverage, prob = p)
  }
  else
  {
    log_p = sapply(
      nvs,
      VGAM::dbetabinom,
      size = coverage,
      prob = p,
      rho = rho
    )
  }

  # Compute P(X > NV) for each NV value
  p_x = cumsum(log_p)

  # Find l_a such that P(X <= l_a) < alpha
  l_a = which(p_x < alpha_level, arr.ind = TRUE) %>% max
  # Find r_a such that P(X > r_a) < 1 - alpha
  r_a = which(p_x > 1 - alpha_level, arr.ind = TRUE) %>% min
  # Adjustments for plots when test fails
  if(is.infinite(l_a)) l_a = 1
  if(is.infinite(r_a)) r_a = coverage

  # Translate NV cutoffs in VAF space
  vafs = nvs / coverage
  l_v = vafs[l_a]
  r_v = vafs[r_a]

  inputs = data.frame(nv = nvs,
                      p = p_x,
                      VAF = vafs)

  return(list(
    model = model,
    density = inputs,
    rho = rho,
    coverage = coverage,
    purity = purity,
    alpha_level = alpha_level,
    nv = c(l_a, r_a),
    vaf = c(l_v, r_v)
  ))
}
