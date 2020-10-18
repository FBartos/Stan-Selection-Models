#include /likelihood_functions.stan
data {
  // the data
  int<lower=0> K;     // number of studies
  vector[K] t;        // observed test-statistic
  vector[K] N;        // number of participants
  
  // support data
  int J;              // number weights cut-offs
  matrix[K,J] crit_t; // critical test-statistic
  
  // values for prior distributions
  real prior_mu1;
  real prior_mu2;
  real prior_tau1;
  real prior_tau2;
  vector[J+1] prior_weight;

  // additional settings
  int test_type;    // type of the test (1 = one sample, 2 = two sample)
}
transformed data{
  vector[K] df;     // degrees of freedom
  if(test_type == 1){
    df = N - 1;
  }else{
    df = N - 2;
  }
}
parameters {
  vector[K] eta;      // unscaled deviation from the mean treatment effect
  real mu;            // mean treatment effect
  real<lower=0> tau;  // deviation of treatment effects
  simplex[J+1] omega; // weight parameters - transformed into beta by cumsum
}
transformed parameters {
  vector[K] theta = mu + tau * eta;   //per-trial treatment effect
  vector[J+1] beta = cumulative_sum(omega);
  vector[K] ncp;
  if(test_type == 1){
    ncp = theta .* sqrt(N);
  }else{
    ncp = theta .* sqrt(N)/2;
  }
}
model {
  // priors
  target += normal_lpdf(mu       | prior_mu1,  prior_mu2);
  target += inv_gamma_lpdf(tau   | prior_tau1, prior_tau2); 
  target += dirichlet_lpdf(omega | prior_weight);
  
  // model
  target += normal_lpdf(eta  | 0, 1); 
  // model
  for(k in 1:K){
    target += w2t_lpdf(t[k] | ncp[k], beta, df[k], crit_t[k,], J);
  }
}

