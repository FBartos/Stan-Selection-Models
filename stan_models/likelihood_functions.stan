functions{
  // denominator likelihood (1)
  real w1tdenom(real x, real xc, real[] param, real[] x_r, int[] x_i){
    
    real ncp    = param[1];
    real beta   = param[2];
    real t      = x_r[1];
    real df     = x_r[2];
    real crit_t = x_r[3];
    real w;
    real lik;
    
    if(fabs(x) < crit_t){
      w = (fabs(x)/crit_t)^beta;
    }else{
      w = 1;
    }
    
    if(ncp > 0){
      lik = exp(student_t_lpdf(x | df, ncp, 1)+log(w));
    }else{
      lik = exp(student_t_lpdf(-x | df, -ncp, 1)+log(w));
    }
    
    return lik;
  }
  // weighted likelihood (1)
  real w1t_lpdf(real t, real ncp, real beta, real df, real crit_t, real[] x_r, int[] x_i, real tolerance){
    
    real log_lik;
    real nominator;
    real denominator;
    real w;
    
    // setting weights
    if(fabs(t) < crit_t){
      w = (fabs(t)/crit_t)^beta;
    }else{
      w = 1;
    }
    
    if(ncp > 0){
      nominator   = student_t_lpdf(t | df, ncp, 1)+log(w);
    }else{
      nominator   = student_t_lpdf(-t | df, -ncp, 1)+log(w);
    }
    
    denominator = log(integrate_1d(w1tdenom, negative_infinity(), positive_infinity(),{ncp, beta}, x_r, x_i, tolerance));
    
    log_lik = nominator-denominator;
    
    return log_lik;
  }
  
  // weighted likelihood (2)
  real w2t_lpdf(real t, real ncp, vector beta, real df, row_vector crit_t, int J){
    
    real log_lik;
    real nominator;
    real denominator;
    real w;
    real denom_part[J+2];
    
    //// setting weights
    // changed to using beta as weights, due to better prior behavior
    if(fabs(t) >= crit_t[J]){
      w = 1;
    }else if(fabs(t) < crit_t[1]){
      w = beta[1];
    }else{
      for(j in 2:J){
        if(fabs(t) < crit_t[j] && fabs(t) >= crit_t[j-1]){
          w = beta[j];
        }
      }
    }
    
    // abusing symetry around zero to mitigate sampling issues for large negative values
    if(ncp > 0){
      
      nominator   = student_t_lpdf(t | df, ncp, 1)+log(w);
      
      // the lower unweighted tail
      denom_part[J+1] = student_t_lcdf(-crit_t[J] | df, ncp, 1);
      // the upper unweighted tail
      denom_part[J+2] = student_t_lccdf(crit_t[J] | df, ncp, 1); 
      
      // the weighted middle
      denom_part[1] = log(exp(student_t_lcdf(crit_t[1] | df, ncp, 1)) - exp(student_t_lcdf(-crit_t[1] | df, ncp, 1)))+log(beta[1]);
      // in an 'onion' slices fasion
      if(J > 1){
        for(jj in 2:J){
          // in an 'onion' slices fasion
          denom_part[jj] = log(
            (exp(student_t_lcdf(-crit_t[jj-1] | df, ncp, 1)) - exp(student_t_lcdf(-crit_t[jj]   | df, ncp, 1))) +
              (exp(student_t_lcdf( crit_t[jj] | df, ncp, 1)) - exp(student_t_lcdf( crit_t[jj-1] | df, ncp, 1)))) +
            log(sum(beta[1:jj]));
        }
      }
  
    }else{
      
      nominator   = student_t_lpdf(-t | df, -ncp, 1)+log(w);
      
      // the lower unweighted tail
      denom_part[J+1] = student_t_lccdf(crit_t[J] | df, -ncp, 1);
      // the upper unweighted tail
      denom_part[J+2] = student_t_lcdf(-crit_t[J] | df, -ncp, 1); 
      
      // the weighted middle
      denom_part[1] = log(exp(student_t_lccdf(-crit_t[1] | df, -ncp, 1)) - exp(student_t_lccdf(crit_t[1] | df, -ncp, 1)))+log(beta[1]);
       if(J > 1){
        for(jj in 2:J){
          // in an 'onion' slices fasion
          denom_part[jj] = log(
            (exp(student_t_lccdf( crit_t[jj-1] | df, -ncp, 1)) - exp(student_t_lccdf(crit_t[jj]    | df, -ncp, 1))) +
              (exp(student_t_lccdf(-crit_t[jj] | df, -ncp, 1)) - exp(student_t_lccdf(-crit_t[jj-1] | df, -ncp, 1)))) +
            log(sum(beta[1:jj]));
        }
      }
      
    }
    
    denominator = log_sum_exp(denom_part);
    
    log_lik = nominator-denominator;
    
    return log_lik;
  }
  
  // denominator likelihood (3 and 4)
  real w3tdenom(real x, real xc, real[] param, real[] x_r, int[] x_i){
    
    real ncp    = param[1];
    real beta   = param[2];
    real t      = x_r[1];
    real df     = x_r[2];
    real power  = x_r[3];
    real pt;
    real lik;
    
    pt  = exp(student_t_lccdf(fabs(x) | df, 0, 1));
    if(ncp > 0){
      lik = exp(student_t_lpdf(x | df, ncp, 1)+(-beta*pt^power));
    }else{
      lik = exp(student_t_lpdf(-x | df, -ncp, 1)+(-beta*pt^power));
    }
    
    return lik;
  }
  // weighted likelihood (3 and 4)
  real w3t_lpdf(real t, real ncp, real beta, real df, real crit_t, real[] x_r, int[] x_i, real tolerance){
    
    real log_lik;
    real nominator;
    real denominator;
    real pt;
    real power  = x_r[3]; //allows to reuse the same functions for weights 3 and 4
    
    // don't try to compute the likelihood for an exteremely unlikely values
    if((fabs(t/ncp) > 5.0 || fabs(t/ncp) < 0.2) && fabs(fabs(t) - fabs(ncp)) > 10.0)return(log(0));
    
    pt = exp(student_t_lccdf(fabs(t) | df, 0, 1));
    if(ncp > 0){
      nominator = student_t_lpdf(t | df, ncp, 1)+(-beta*pt^power);
    }else{
      nominator = student_t_lpdf(-t | df, -ncp, 1)+(-beta*pt^power);
    }
    denominator = log(integrate_1d(w3tdenom, negative_infinity(), positive_infinity(),{ncp, beta}, x_r, x_i, tolerance));
    
    log_lik = nominator-denominator;
    
    return log_lik;
  }
}
