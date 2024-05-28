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
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] smax_dev; //year-to-year deviations in Smax

}

transformed parameters{
  vector[L] Smax; //Smax in each year
  vector[L] b; //rate capacity in each year
  
  Smax[1] = 1/b0;
  for(t in 2:L){
    Smax[t] = Smax[t-1] + smax_dev[t-1]*sigma_b;
  }
 b=1/Smax;
}  

model{
  //priors
  log_a ~ normal(1.5,2.5); //productivity
  b0 ~ lognormal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter
  smax_dev ~ std_normal();

  //variance terms
  sigma ~ normal(0.5,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,pSmax_sig); //half normal on variance (lower limit of zero)
  
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma);
}
generated quantities{
     real Umsy;
     vector[L] Smsy;
         
    for(l in 1:L){ Smsy[l] = (1-lambert_w0(exp(1-log_a)))/b[l];
    }
    Umsy = 1-lambert_w0(exp(1-log_a));
}
 
