compute_density = function(NV, DP, prob, model, rho){
  if ((model %>% tolower()) == 'binomial') {
    density = dbinom(x = 1:DP, size =  DP, prob =  prob)
  }
  else {
    density = VGAM::dbetabinom(x = 1:DP, size = DP, prob = prob, rho = rho) 
    }
}
