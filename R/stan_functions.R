#' sr_mod function
#'
#' This function generates a stock-recruitment (S-R) model for rstan based on user inputs.
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'tv' model, or a regime shift model 'regime'
#' @param ac TRUE or FALSE statement to include autocorrelated residuals. Only compatible with static model
#' @param par For time-varying or regime S-R models, what parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param loglik TRUE or FALSE statement that dictates whether model is being used for out-of-sample log-likelihood estimation
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model()
#' @export
#' @examples
#' m2=sr_mod(type='static',ac = TRUE,par='n',loglik=T)

sr_mod<- function(type=c('static','tv','regime'),ac=FALSE,par=c('n','a','b','both'),loglik=FALSE){
  if(type=='static'&ac==F){
    if(loglik==FALSE){
      m="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<upper = 0> log_b; // rate capacity - fixed in this
    
     //variance components  
      real<lower = 0> sigma_e;
    
    }
    transformed parameters{
    	real b;
    	
    	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
      //priors
      log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
      sigma_e ~ gamma(2,3);
      
      R_S ~ normal(log_a - S*b, sigma_e);
    }
    "
    }
    if(loglik==TRUE){
      m ="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      real y_oos; //log(recruits per spawner)
      real x_oos; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<upper = 0> log_b; // rate capacity - fixed in this
    
     //variance components  
      real<lower = 0> sigma_e;
    
    }
    transformed parameters{
    	real b;
    	
    	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
      //priors
      log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
      sigma_e ~ gamma(2,3);
      
      R_S ~ normal(log_a - S*b, sigma_e);
    }
    generated quantities{
     real log_lik_oos;
    	log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma_e);
    }
    "}
  }
  if(type=='static'&ac==T){
    if(loglik==FALSE){
      m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> TT;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;

}
transformed parameters{
real b;
vector[N] mu;
vector[N] epsilon; //residuals

b = exp(log_b);
mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (phi^(ii[t]-ii[t-1])*epsilon[t-1]); //phi raised the power of the number of time-steps between successive productivity estimates
  }
}
model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ gamma(2,3);

  R_S ~ normal(mu, sigma_e);
  
}
    "
    }
if(loglik==TRUE){
  m ="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T

 }
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;

}
transformed parameters{
real b;
vector[N] mu;
vector[N] epsilon; //residuals

b = exp(log_b);
mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (phi^(ii[t]-ii[t-1])*epsilon[t-1]);
  }
}
model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ gamma(2,3);

  R_S ~ normal(mu, sigma_e);
  
}
generated quantities{
  real ep_3b;
  real ep_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
   
  ep_3b = (epsilon[N]+epsilon[N-1]+epsilon[N-2])/3;
  ep_5b = (epsilon[N]+epsilon[N-1]+epsilon[N-2]+epsilon[N-3]+epsilon[N-4])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b+phi*epsilon[N], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b+phi*ep_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b+phi*ep_5b, sigma_e);
 } 
    "}
  }
if(type=='tv'&par=='a'){
  if(loglik==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
  
}  
model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
   
 
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma_e); 
  
} "
  }
if(loglik==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
}  
model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
   
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma_e);
}
  generated quantities{
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b, sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma_e);
 }
    "}
}  
if(type=='tv'&par=='b'){
  if(loglik==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma_e);
}
 "
  }
if(loglik==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a-S[n]*b[ii[n]], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma_e);
 }
 "}
}
if(type=='tv'&par=='both'){
  if(loglik==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}"
  }
if(loglik==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
   
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_e);
 }
 "}
}
if(type=='regime'&par=='a'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator1[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator2[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator1[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real log_a_1b;
real sigma_1b;
real log_a_3b;
real sigma_3b;
real log_a_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] log_a_1bw_k;
vector[K] sigma_1bw_k;
vector[K] log_a_3bw_k;
vector[K] sigma_3bw_k;
vector[K] log_a_5bw_k;
vector[K] sigma_5bw_k;

real log_a_1bw;
real sigma_1bw;
real log_a_3bw;
real sigma_3bw;
real log_a_5bw;
real sigma_5bw;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator2[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
sigma_1b = sigma[zstar[N]]; //sigma error
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

//slope weighted by probability of each regime 
for(k in 1:K) log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime

log_a_1bw=sum(log_a_1bw_k); //weighted productivity
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
log_a_3bw=sum(log_a_3bw_k); //weighted productivity
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
log_a_5bw=sum(log_a_5bw_k); //weighted productivity
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a_1bw - x_oos*b, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a_3bw - x_oos*b, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a_5bw - x_oos*b, sigma_5bw);
}
 "}
}
if(type=='regime'&par=='b'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K]<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,2.5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];
{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K] log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,2.5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real b_1b;
real sigma_1b;
real b_3b;
real sigma_3b;
real b_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] b_1bw_k;
vector[K] sigma_1bw_k;
vector[K] b_3bw_k;
vector[K] sigma_3bw_k;
vector[K] b_5bw_k;
vector[K] sigma_5bw_k;

real b_1bw;
real sigma_1bw;
real b_3bw;
real sigma_3bw;
real b_5bw;
real sigma_5bw;


{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

b_1b = b[zstar[N]]; //slope based on most probable state in sample N
sigma_1b = sigma[zstar[N]]; //sigma error of most probable state in sample N
b_3b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
b_5b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]]+b[zstar[N-3]]+b[zstar[N-4]])/5; //intercept
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

for(k in 1:K){
 b_1bw_k[k]=gamma[N,k]*b[k]; //prob of each regime x productivity for each regime
 sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime
 b_3bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k])/3; //prob of each regime x productivity for each regime
 sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime
 b_5bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k]+gamma[N-3,k]*b[k]+gamma[N-4,k]*b[k])/5; //prob of each regime x productivity for each regime
sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime
}

b_1bw=sum(b_1bw_k); //weighted capacity
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
b_3bw=sum(b_3bw_k); //weighted productivity
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
b_5bw=sum(b_5bw_k); //weighted productivity
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b_1b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a - x_oos*b_1bw, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a - x_oos*b_3bw, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a - x_oos*b_5bw, sigma_5bw);
}"}
}
if(type=='regime'&par=='both'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K]<lower = 0>  log_a; // regime max. productivity
ordered[K]<upper = 0> log_b; // regime rate capacity 
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];
{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K]  log_a; // regime max. productivity
ordered[K] log_b; // regime rate capacity 
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real log_a_1b;
real log_a_3b;
real log_a_5b;
real b_1b;
real sigma_1b;
real b_3b;
real sigma_3b;
real b_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] log_a_1bw_k;
vector[K] sigma_1bw_k;
vector[K] log_a_3bw_k;
vector[K] sigma_3bw_k;
vector[K] log_a_5bw_k;
vector[K] sigma_5bw_k;
vector[K] b_1bw_k;
vector[K] b_3bw_k;
vector[K] b_5bw_k;

real b_1bw;
real b_3bw;
real b_5bw;
real log_a_1bw;
real sigma_1bw;
real log_a_3bw;
real sigma_3bw;
real log_a_5bw;
real sigma_5bw;


{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
b_1b = b[zstar[N]]; //slope based on most probable state in sample N
sigma_1b = sigma[zstar[N]]; //sigma error
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept 3-y back
b_3b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
b_5b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]]+b[zstar[N-3]]+b[zstar[N-4]])/5; 
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

//slope weighted by probability of each regime
//slope weighted by probability of each regime 
for(k in 1:K){
log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime
b_1bw_k[k]=gamma[N,k]*b[k]; //prob of each regime x productivity for each regime
sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime

log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime
b_3bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k])/3; //prob of each regime x productivity for each regime
sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime

log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime
b_5bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k]+gamma[N-3,k]*b[k]+gamma[N-4,k]*b[k])/5; //prob of each regime x productivity for each regime
sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime
}
log_a_1bw=sum(log_a_1bw_k); //weighted productivity
b_1bw=sum(b_1bw_k); //weighted capacity - 1 year previous
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
log_a_3bw=sum(log_a_3bw_k); //weighted productivity
b_3bw=sum(b_3bw_k); //weighted capacity - 3 year previous average
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
log_a_5bw=sum(log_a_5bw_k); //weighted productivity
b_5bw=sum(b_5bw_k); //weighted capacity - 5 year previous average
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b_1b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a_1bw - x_oos*b_1bw, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a_3bw - x_oos*b_1bw, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a_5bw - x_oos*b_5bw, sigma_5bw);

}

"}
}
m2=rstan::stan_model(model_code = m)
return(m2)  
}
