#include /likelihood_functions.stan
data {
  // the data
  int<lower=0> K;     // number of studies
  vector[K] t;        // observed test-statistic
  vector[K] N;        // number of participants

  // support data
  vector[K] crit_t;   // critical test-statistic
  
  // values for prior distributions
  real prior_mu1;
  real prior_mu2;
  real prior_weight1;
  real prior_weight2;
  
  int test_type;    // type of the test (1 = one sample, 2 = two sample)
  int weight_type;  // type of the likelihood weight
  real tolerance;   // integration tolerance settings
}
transformed data{
  int x_i[0];       // empty object just for weighted likelihood integration
  vector[K] df;
  if(test_type == 1){
    df = N - 1;
  }else{
    df = N - 2;
  }
}
parameters{
  real mu;            // mean treatment effect
  real<lower=0> beta; // weight coefficient
}
transformed parameters{
  vector[K] ncp = mu*sqrt(N)/2;
}
model{
  // priors
  target += normal_lpdf(mu      | prior_mu1, prior_mu2);
  target += lognormal_lcdf(beta | log(prior_weight1), log(prior_weight2));
  
  // model
  if(weight_type == 1){
    for(k in 1:K){
      target += w1t_lpdf(t[k] | ncp[k], beta, df[k], crit_t[k], {t[k], df[k], crit_t[k]}, x_i, tolerance);
    }
  }else if(weight_type == 3){
    for(k in 1:K){
      target += w3t_lpdf(t[k] | ncp[k], beta, df[k], crit_t[k], {t[k], df[k], 2}, x_i, tolerance);
    }
  }else if(weight_type == 4){
    for(k in 1:K){
      target += w3t_lpdf(t[k] | ncp[k], beta, df[k], crit_t[k], {t[k], df[k], 1}, x_i, tolerance);
    }
  }
}
