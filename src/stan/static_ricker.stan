data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
transformed data{
real logbeta_pr;
real logbeta_pr_sig;

logbeta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr=log(1/pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction
}
parameters {
  real log_a;// initial productivity (on log scale)
  real<lower = 0> b; // rate capacity
    
//variance components  
  real<lower = 0> sigma;
    
}
model{
  //priors
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  b ~ lognormal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter - wide prior
   
  //variance terms
  sigma ~ normal(0.5,1); //half-normal prior - expectation around 0.5

   R_S ~ normal(log_a - S*b, sigma);
}
generated quantities{
 real Smax;
 real Umsy;
 real Smsy;
 
Smax = 1/b;
Umsy = 1-lambert_w0(exp(1-log_a));
Smsy = (1-lambert_w0(exp(1-log_a)))/b;
}
    
