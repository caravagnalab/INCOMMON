functions {
  real lambda(int k, real x, real purity) {
    return 2 * x * (1 - purity) + k * x * purity;
  }

  real expected_vaf(int m, int k, real purity) {
    return m * purity / (2 * (1 - purity) + k * purity);
  }

  vector alpha_beta_pi(real purity_mean, real purity_error) {
    real purity_variance = square(purity_error);
    real alpha_pi = purity_mean * ((purity_mean * (1 - purity_mean) / purity_variance) - 1);
    real beta_pi = (1 - purity_mean) * ((purity_mean * (1 - purity_mean) / purity_variance) - 1);
    vector[2] alpha_beta;
    alpha_beta[1] = alpha_pi;
    alpha_beta[2] = beta_pi;
    return alpha_beta;
  }
}

data {
  int<lower=0> M;          // number of mutations
  array[M] int<lower=0> N; // observed counts
  array[M] int<lower=0> n; // variant read counts for each mutation
  int<lower=1> k_max;
  real<lower=0, upper=1> purity_mean;
  real<lower=0, upper=1> purity_error;
  real<lower=0> alpha_x;  // Gamma prior parameter alpha for x
  real<lower=0> beta_x;   // Gamma prior parameter beta for x
  array[M] vector[k_max*(k_max+1)/2] alpha_k_m;  // Dirichlet priors for each mutation
}

parameters {
  array[M] simplex[k_max*(k_max+1)/2] psi;  // Mixing proportions for each mutation
  real<lower=0> x;    // Per-copy count rate
  real<lower=0,upper=1> purity; // Proportion of tumour cells (purity)
}

model {
  // Prior sampling
  vector[2] ab_pi = alpha_beta_pi(purity_mean, purity_error); // Compute shape parameters for purity prior
  purity ~ beta(ab_pi[1], ab_pi[2]);  // Prior sampling of the purity
  x ~ gamma(alpha_x, beta_x);  // Prior sampling of the per-count rate

  for (i in 1:M) {
    psi[i] ~ dirichlet(alpha_k_m[i]);  // Prior sampling of each mutation's mixing proportions
  }

  // Likelihood sampling
  for (i in 1:M) {
    vector[k_max*(k_max+1)/2] log_likelihood;
    int idx = 1;
    for (k in 1:k_max) {
      for (m in 1:k) {
        log_likelihood[idx] = binomial_lpmf(n[i] | N[i], expected_vaf(m, k, purity))
          + poisson_lpmf(N[i] | lambda(k, x, purity));
        idx += 1;
      }
    }
    target += log_sum_exp(log_likelihood + log(psi[i]));
  }
}
