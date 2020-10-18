data {
  // the data
  int<lower=0> K;     // number of studies
  vector[K] t;        // observed test-statistic
  vector[K] N;        // number of participants
  
  // values for prior distributions
  real prior_mu1;
  real prior_mu2;

  // additional settings
  int test_type;      // type of the test (1 = one sample, 2 = two sample)
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
  real mu;            // mean treatment effect
}
transformed parameters{
  vector[K] ncp = mu*sqrt(N)/2;
}
model {
  // priors
  target += normal_lpdf(mu   | prior_mu1, prior_mu2);
  
  // model
  target += student_t_lpdf(t | df, ncp, 1); 
}
