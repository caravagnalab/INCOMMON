compute_likelihood = function(NV, DP, prob, rho) {
  density = VGAM::dbetabinom(
    x = 1:DP,
    size = DP,
    prob = prob,
    rho = rho
  )
}
