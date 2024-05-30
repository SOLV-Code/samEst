data{
  int<lower=1> N;//number of annual samples
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
transformed data{
real logbeta_pr;
real logbeta_pr_sig;

logbeta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr=log(1/pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction
}
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity (on log scale)

  //variance components  
  real<lower=0> sigma; //observation error at t=1
  real<lower = -1, upper = 1> rho; //autocorrelation parameter

}
transformed parameters{
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  real<lower=0> sigma_AR; //observation error - corrected for autocorrelation
  real<lower=0> b=exp(log_b); //per capita density dependence
  
  mu = log_a-b*S; //inital expectation

  epsilon[1] = R_S[1] - mu[1]; //first observation without observed autocorrelation
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]); //residual productivity
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]); //rho raised the power of the number of time-steps between successive productivity estimates
  }
  sigma_AR = sigma*sqrt(1-rho^2); //correct observation error for autocorrelation
}
model{
  //priors
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(logbeta_pr,logbeta_pr_sig); //per capita capacity parameter - wide prior
     
  //variance terms
  sigma ~ normal(0.5,1); //half-normal prior - expectation at 0.5
  //autocorrelation term
  rho ~ uniform(-1,1);
  
  R_S[1] ~ normal(mu[1], sigma);
  for(t in 2:N) R_S[t] ~ normal(mu[t], sigma_AR);
  
}
generated quantities{
  real Smax;
  real Umsy;
  real Smsy;
    
  Smax = 1/b;
  Umsy = 1-lambert_w0(exp(1-log_a));
  Smsy = (1-lambert_w0(exp(1-log_a)))/b;
}
    
