data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
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
  real log_a;// initial productivity (on log scale) - fixed in this
  real log_b0; // inital log rate capacity

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in Smax

}

transformed parameters{
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  vector[L] b; //rate capacity in each year
  
  log_b[1] = 1/log_b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
    }
b=exp(log_b);
 mu=log_a-b[ii]*S;
  epsilon=R_S-mu;
}  

model{
  //priors
  log_a ~ normal(1.5,2.5); //productivity
  log_b0 ~ lognormal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter
  b_dev ~ std_normal();

  //variance terms
  sigma ~ normal(0.5,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,1); //half normal on variance (lower limit of zero)
  
  R_S ~ normal(mu, sigma);
}
generated quantities{
     real Umsy;
     vector[L] Smsy;
	vector[L] Smax;
         
    for(l in 1:L){ 
    Smsy[l] = (1-lambert_w0(exp(1-log_a)))/b[l];
    Smax[l] = 1/b[l];
}
    Umsy = 1-lambert_w0(exp(1-log_a));
}
 
