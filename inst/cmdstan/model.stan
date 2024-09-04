functions{

  real lambda(int k, real x, real purity){
    return 2 * x * (1 - purity) + k * x * purity;
  }
  
  real expected_vaf(int m, int k, real purity){
    return m * purity / (2 * (1 - purity) + k * purity);
  }
  
  vector alpha_beta_pi(real purity_mean, real purity_error) {
    
    real purity_variance = square(purity_error);
    
    real alpha_pi = purity_mean * ((purity_mean * (1 - purity_mean) / purity_variance) - 1);
    real beta_pi = (1 - purity_mean) * ((purity_mean * (1 - purity_mean) / purity_variance) - 1);
    
    // Return a vector containing alpha and beta
    vector[2] alpha_beta;
    alpha_beta[1] = alpha_pi;
    alpha_beta[2] = beta_pi;

    return alpha_beta;
  }
  
}

data {
  int<lower=0> M;         // number of mutations
  array[M] int<lower=0> N;     // observed counts
  array[M] int<lower=0> n;    // variant read counts for each mutation
  int<lower=1> k_max;
  real<lower=0, upper=1> purity_mean;
  real<lower=0, upper=1> purity_error;
  real<lower=0> alpha_x;  // Gamma prior parameter alpha for x
  real<lower=0> beta_x;   // Gamma prior parameter beta for x
  // array[M] simplex[K_max] ploidy_prior; 
  array[M] simplex[k_max*k_max] k_m_prior; 
}

parameters {
  real<lower=0> x;    // rate for Poisson distribution (x)
  real<lower=0,upper=1> purity; // Proportion of tumour cells (purity)
  // Probability of each multiplicity class for each mutation
  array[M] simplex[3] class_probs;  // Mixing proportions for each mutation
}

model {

  vector[2] ab_pi = alpha_beta_pi(purity_mean, purity_error);
  purity ~ beta(ab_pi[1], ab_pi[2]); // Beta prior for pi with the new parameters
  x ~ gamma(alpha_x, beta_x);  // Gamma prior for x

  for (i in 1:M) {
    
    // Compute the joint log likelihood of n | x, class = (m, k), purity and N | x, k, purity
    array[k_max, 3] real log_lik;
    
    // Class m = 1 (k > 1)
    log_lik[1, 1] = negative_infinity();
    for(k in 2:k_max){
      // log_lik[k, 1] = binomial_lpmf(n[i] | N[i], expected_vaf(1, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(ploidy_prior[i, k]);
      log_lik[k, 1] = binomial_lpmf(n[i] | N[i], expected_vaf(1, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(k_m_prior[i, (k-1)*k_max+1]);
    }
    
    // Class m = k
    for(k in 1:k_max){
      // log_lik[k, 2] = binomial_lpmf(n[i] | N[i], expected_vaf(k, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(ploidy_prior[i, k]);
      log_lik[k, 2] = binomial_lpmf(n[i] | N[i], expected_vaf(k, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(k_m_prior[i, (k-1)*k_max+k]);
    }
    
     // Class 1 < m < k (k > 2)
    log_lik[1, 3] = negative_infinity();
    log_lik[2, 3] = negative_infinity();
    for(k in 3:k_max){
      array[k - 2] real lps;
      for(m in 2:(k-1)){
        // lps[m-1] = binomial_lpmf(n[i] | N[i], expected_vaf(m, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(ploidy_prior[i, k]);
        lps[m-1] = binomial_lpmf(n[i] | N[i], expected_vaf(m, k, purity)) + poisson_lpmf(N[i] | lambda(k, x, purity)) + log(k_m_prior[i, (k-1)*k_max+m]);
      }
      log_lik[k, 3] = log_sum_exp(lps);
    }
    
    // Compute joint log likelihood of classes (m, k)
    array[3] real log_posterior_class;
    log_posterior_class[1] = log_sum_exp(log_lik[, 1]) + log(class_probs[i, 1]);
    log_posterior_class[2] = log_sum_exp(log_lik[, 2]) + log(class_probs[i, 2]);
    log_posterior_class[3] = log_sum_exp(log_lik[, 3]) + log(class_probs[i, 3]);
    
    target += log_sum_exp(log_posterior_class);
    
  }
}
