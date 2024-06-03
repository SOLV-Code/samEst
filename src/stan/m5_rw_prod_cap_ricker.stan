data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  matrix[1,N] S; //spawners in time T
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
  real log_a0;// initial productivity (on log scale)
  real log_b0; // initial rate capacity (log)

 //variance components  
  real<lower = 0> sigma;
 real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //z-score annual deviations in log_a
  vector[L-1] b_dev; //z-score annual deviations in Smax

}

transformed parameters{
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //log rate capacity in each year
  vector[L] b; //rate capacity in each year
  
  log_a[1] = log_a0;
  log_b[1] = 1/log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  }
 b=exp(log_b);
 mu=log_a[ii]-S*b[ii];
 epsilon=R_S-mu;
}  

model{
  //priors
  log_a0 ~ normal(1.5,5); //productivity
  log_b0 ~ normal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();

  //variance terms
  sigma ~ normal(0.5,1); //half normal on variance
  sigma_a ~ normal(0,1); //half normal on variance 
  sigma_b ~ normal(0,1); //informative prior for Smax so it knows the parameter scale
  
 R_S ~ normal(mu, sigma);
}
generated quantities{
     vector[L] Umsy;
     vector[L] Smsy;
     vector[L] Smax;

    for(l in 1:L){ 
Smsy[l] = (1-lambert_w0(exp(1-log_a[l])))/b[l];
Smax[l] = 1/b[l];
  Umsy[l] = 1-lambert_w0(exp(1-log_a[l]));
    }
  
}
 
