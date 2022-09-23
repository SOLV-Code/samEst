#include <TMB.hpp>
#include <iostream>


template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> a, int do_log)
{
  Type phi = a.sum();
  int n = x.size();
  Type ll = lgamma(phi);
  for(int i = 0; i < n; i++) ll +=  -lgamma(a(i)) + (a(i) - 1.0) * log(x(i));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type ddirmultinom (vector<Type> obs, vector<Type> p, Type phi, int do_log) {
  int dim = obs.size();
  Type N = obs.sum();
  // e.g. Eq. 4 Thorson et al. 2017 Fish. Res.; phi here = Beta in paper
  // or Eq. B.2 in https://www.sciencedirect.com/science/article/pii/S0165783621000953#sec0190
  // this is the 'saturating' version
  Type ll = lgamma(N + Type(1.0)) + lgamma(phi) - lgamma(N + phi); // eq. outside of summations
  for (int a = 0; a < dim; a++) {
    ll += -lgamma(obs(a) + Type(1.0)) +
      lgamma(obs(a) + phi * (p(a) + Type(1.0e-15))); // 1e-15 for robustness to 0s
      lgamma(phi * (p(a) + Type(1.0e-15)));
  }
  if (do_log) return ll;
  else return exp(ll);
}


template<class Type>
vector<Type> segment_1(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type> pi1,vector<Type> alpha,vector<Type> beta,Type sigma,int t){
  
  int k_regime = beta.size();
  Type small = pow(10,-300);
  vector<Type> sr = log(pi1 + small);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = alpha(j) - beta(j)*st(0);
    sr(j) += dnorm(yt(0), f_now, sigma,true);
  }
  
  for(int i = 1;i <= t;++i){
    vector<Type> sr_new = sr;    
    
    for(int j = 0;j < k_regime;++j){
      sr_new(j) = sr(0) +qij(0,j);
    
      for(int jj = 1;jj < k_regime;++jj){
        Type temp = sr(jj) +qij(jj,j);
        sr_new(j) = logspace_add(sr_new(j),temp);
      }
    }
    sr = sr_new;
 
    for(int j = 0;j < k_regime;++j){
      Type f_now = alpha(j) - beta(j)*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma,true);
    }
  }
  return sr;
}


template<class Type>
Type segment_2(vector<Type> yt, vector<Type> st, matrix<Type> qij,vector<Type>
pi1,vector<Type> alpha,vector<Type> beta,Type sigma,int rt,int t){
  
  int k_regime = beta.size();
  int n = yt.size();
  vector<Type> sr = qij.row(rt);
  
  for(int j = 0;j < k_regime;++j){
    Type f_now = alpha(j) - beta(j)*st(t+1);
    sr(j) += dnorm(yt(t+1), f_now, sigma,true);
  }
 
  for(int i = t+2;i < n;++i){
    vector<Type> sr_new = sr;
    for(int j = 0;j < k_regime;++j){
      sr_new(j) = sr(0) +qij(0,j);
      for(int jj = 1;jj < k_regime;++jj){
        Type temp = sr(jj) +qij(jj,j);
        sr_new(j) = logspace_add(sr_new(j),temp);
      }
    }
    sr = sr_new;
    for(int j = 0;j < k_regime;++j){
      Type f_now = alpha(j) - beta(j)*st(i);
      sr(j) += dnorm(yt(i), f_now, sigma,true);
    }
  }
  Type seg2 = sr(0);
  for(int j = 1;j < k_regime;++j){
    seg2 = logspace_add(seg2,sr(j));
  }
  return seg2;
}

//main //////////////////////////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(yt);
  DATA_VECTOR(st);
  DATA_SCALAR(alpha_u); //upper bound for b
  DATA_SCALAR(alpha_l); //lower bound for b
  DATA_SCALAR(beta_u);  //upper bound for a
  DATA_SCALAR(sigma_u); //upper bound for sigma
  DATA_VECTOR(alpha_dirichlet); //prior inputs for dirichlet 

  PARAMETER_VECTOR(lalpha);
  PARAMETER_VECTOR(lbeta);
  PARAMETER(lsigma);
  PARAMETER_VECTOR(pi1_tran); // initial state probabilities
  PARAMETER_MATRIX(qij_tran); // transition probabilities

  int k_regime = lbeta.size();
  
  vector<Type> beta = beta_u/(1+exp(-lbeta));// when lbeta is negative infinity, beta=0; when lbeta is positive infinity, beta=beta_u
  vector<Type> alpha(k_regime);
  
  alpha(0) = (alpha_u-alpha_l)/(1+exp(-lalpha(0)))+alpha_l;
  
  for(int i = 1;i < k_regime;++i){
    alpha(i) = alpha(i-1) + (alpha_u-alpha(i-1))/(1+exp(-lalpha(i)));
  } // alpha(1) from alpha(0) to alpha_u

  Type sigma = sigma_u/(1+exp(-lsigma));
  vector<Type> pi1(k_regime); 
 
  for(int i = 0;i < k_regime-1;++i){
    pi1(i) = exp(pi1_tran(i));    
  }
  pi1(k_regime-1) = 1;
  pi1 = pi1/(pi1.sum());
  
  Type small = pow(10,-300);
  matrix<Type> qij(k_regime,k_regime);
  for(int i = 0;i < k_regime;++i){
    for(int j = 0;j < k_regime-1;++j){
      qij(i,j) = exp(qij_tran(i,j));
    }
    qij(i,k_regime-1) = 1;
    vector<Type> qij_row = qij.row(i);
    Type row_sum = qij_row.sum();
    for(int j = 0;j < k_regime;++j){
      qij(i,j) = qij(i,j)/row_sum;
      qij(i,j) = log(qij(i,j)+small);
    }
  } 
  int n = yt.size();
  vector<Type> sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1);
  Type nll = sr(0);
  
  for(int j = 1;j < k_regime;++j){
    nll = logspace_add(nll,sr(j));
  }
  nll = -nll;

  
  
// predict r_t ///////////////////////////////////
  matrix<Type> r_pred(k_regime,n);
  for(int i = 0;i < n-1;++i){
    sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,i);
    for(int j = 0;j < k_regime;++j){
      Type tempt = sr(j) + segment_2(yt, st, qij,pi1,alpha,beta,sigma,j,i) + nll;
      r_pred(j,i) = exp(tempt);
    }
  }
  sr = segment_1(yt, st, qij,pi1,alpha,beta,sigma,n-1);
  for(int j = 0;j < k_regime;++j){
    Type tempt = sr(j) + nll;
    r_pred(j,n-1) = exp(tempt);
  }
 
 qij = exp(qij.array());


  //priors
  Type pnll = Type(0.0);
  vector<Type> pi_prior(k_regime);
 
  pnll -= dgamma(sigma,Type(2.0),Type(1.0)/Type(3.0),true);
   
  for(int j = 0;j < k_regime;++j){
    pi_prior(j) = Type(1.0);
    Type logbeta = log(beta(j));
    pnll -= dnorm(alpha(j),Type(0.0),Type(2.5),true);
    pnll -= dnorm(logbeta,Type(-12.0),Type(3.0),true);
    
    vector<Type> qijtmp = qij.row(j);
    pnll -= ddirichlet(qijtmp,alpha_dirichlet,true);   
  
  //
  }
  pnll -=ddirichlet(pi1,pi_prior,true);



  Type ans= nll + pnll;

REPORT(beta);
REPORT(alpha);
REPORT(sigma);
REPORT(pi1);
REPORT(qij);
REPORT(r_pred); 
REPORT(nll);
REPORT(pnll);  


 


ADREPORT(alpha);
ADREPORT(beta);
ADREPORT(sigma);
ADREPORT(pi1);
ADREPORT(qij);

return ans;

}
