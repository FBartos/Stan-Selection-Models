data {
  // the data
  int<lower=0> K;     // number of studies
  vector[K] t;        // observed test-statistic
  vector[K] N;        // number of participants
  
  // values for prior distributions
  real prior_mu1;
  real prior_mu2;
  real prior_tau1;
  real prior_tau2;
  
  // additional settings
  int test_type;    // type of the test (1 = one sample, 2 = two sample)
}
transformed data{
  vector[K] df;
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
}
transformed parameters {
  vector[K] theta = mu + tau * eta;   //per-trial treatment effect
  vector[K] ncp;
  if(test_type == 1){
    ncp = theta .* sqrt(N);
  }else{
    ncp = theta .* sqrt(N)/2;
  }
}
model {
  // priors
  target += normal_lpdf(mu     | prior_mu1,  prior_mu2);
  target += inv_gamma_lpdf(tau | prior_tau1, prior_tau2); 
  
  // model
  target += normal_lpdf(eta  | 0, 1); 
  target += student_t_lpdf(t | df, ncp, 1);
}
