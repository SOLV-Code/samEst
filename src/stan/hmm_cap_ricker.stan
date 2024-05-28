functions {
  vector normalize(vector x) {
  return x / sum(x);
}
}
data {
   int<lower=1> N;//number of annual samples
  int<lower=1> L;//time-series length
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet
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
  // Discrete state model
  simplex[K] A[K]; // transition probabilities
 simplex[K] pi1; // initial state probabilities

  // A[i][j] = p(z_t = j | z_{t-1} = i)
  // Continuous observation model
  real log_a; // max. productivity
  ordered<lower=0>[K] b; // rate capacity - fixed in this
  real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
  vector[K] logalpha[N]; //cumulative likelihood (forward)

{ // Forward algorithm log p(z_t = j | y_{1:t})
  real accumulator1[K];

  for(k in 1:K) logalpha[1,k] = log(pi1[k]) + normal_lpdf(R_S[1] |log_a - b[k]*S[1], sigma);

//forward pass:
  for (t in 2:N) {
  for (j in 1:K) { // j = current (t)
  for (i in 1:K) { // i = previous (t-1)
  // Murphy (2012) p. 609 eq. 17.48
  // belief state + transition prob + local evidence at t
    accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma);
  }
  logalpha[t, j] = log_sum_exp(accumulator1);
  }
  }
  }
}
model{
  log_a ~ normal(1.5,2.5);
  for(k in 1:K) b[k] ~ lognormal(logbeta_pr,logbeta_pr_sig); //capacity

  sigma ~ normal(0.5,1); //half normal on variance (lower limit of zero)
  pi1~ dirichlet(rep_vector(1,K));

  for(k in 1:K){
  A[k,] ~ dirichlet(alpha_dirichlet[k,]);
  }
  
  target += log_sum_exp(logalpha[N]);
}
generated quantities {
  vector[N] log_lik;
  int<lower=1, upper=K> zstar[N];
  real logp_zstar;
  vector[K] alpha[N];
  vector[K] logbeta[N];
  vector[K] loggamma[N];
  vector[K] beta[N];
  vector[K] gamma[N];
  
  vector[K] Smax;
  real Umsy;
  vector[K] Smsy;
  
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

  accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a - b[i]*S[t], sigma);
  }
  logbeta[t-1, j] = log_sum_exp(accumulator2);
  }
  }
  for (t in 1:N)
  beta[t] = softmax(logbeta[t]);
  } // Backward

  { // Forward-backward algorithm log p(z_t = j | y_{1:N})
  for(t in 1:N) {
  loggamma[t] = alpha[t] .* beta[t];
  }
  for(t in 1:N)
  gamma[t] = normalize(loggamma[t]);
  } // Forward-backward
  
  { // Viterbi algorithm
  int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
  real delta[N, K]; // max prob for the sequence up to t
  // that ends with an emission from state k
  for (j in 1:K)
  delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma);
  for (t in 2:N) {
    for (j in 1:K) { // j = current (t)
      delta[t, j] = negative_infinity();
      for (i in 1:K) { // i = previous (t-1)
        real logp;
        logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma);
         
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


for(n in 1:N)log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b[zstar[n]], sigma);
   

Umsy = 1-lambert_w0(exp(1-log_a));
for(k in 1:K){
Smax[k] = 1/b[k];
Smsy[k] = (1-lambert_w0(exp(1-log_a)))/b[k];
}

}

