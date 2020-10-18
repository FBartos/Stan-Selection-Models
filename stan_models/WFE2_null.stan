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
  simplex[J+1] omega; // weight parameters - transformed into beta by cumsum
}
transformed parameters {
  vector[J+1] beta = cumulative_sum(omega);
}
model {
  // priors
  target += dirichlet_lpdf(omega | prior_weight);

  // model
  for(k in 1:K){
    target += w2t_lpdf(t[k] | 0, beta, df[k], crit_t[k,], J);
  }
}
