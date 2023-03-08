#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
// Scalar version
template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  

template <class Type>
Type minus_one_to_one(Type x)
{
  return Type(2) * invlogit(x) - Type(1);
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_logRS);   // observed log recruitment
  DATA_INTEGER(priors_flag); //flag indicating wether or not priors should be used
  DATA_INTEGER(stan_flag); //flag indicating wether or not 

  //DATA_SCALAR(sig_p_mean);
  DATA_SCALAR(sig_p_sd); //sd for sigma prior
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  PARAMETER(rho);
  //PARAMETER_VECTOR(delta);
  
  int timeSteps=obs_logRS.size();

  Type rhoo = minus_one_to_one(rho);

  
  Type beta = exp(logbeta);
  Type sigobs = exp(logsigobs);
  Type Smax  = Type(1.0)/beta;
  

  Type sigAR  = sigobs*sqrt(1-pow(rhoo,2));

  
  //priors - based on evaluation done with the prior predictive check
  //Type ans = Type(0);
  Type nll= Type(0);
  Type pnll = Type(0.0);

  if(priors_flag == 1){
    //ans -=dnorm(alpha,Type(0.0),Type(2.5),true);
    pnll -=dgamma(alpha,Type(3.0),Type(1.0),true);
    pnll -=dnorm(logbeta,Type(-12.0),Type(3.0),true); 
    //ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
    //pnll -= dgamma(sigobs,Type(2.0),Type(1.0)/Type(3.0),true);
    pnll -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    if(stan_flag) pnll -= logsigobs;
  }
  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps) ;
 
  pred_logRS(0) = alpha - beta * obs_S(0) ;
  pred_logR(0) = pred_logRS(0) + log(obs_S(0));
  residuals(0) = obs_logRS(0) - pred_logRS(0);
    
  nll+= -dnorm(obs_logRS(0),pred_logRS(0),sigobs,true);

  for(int i=1;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      pred_logRS(i) = alpha - beta * obs_S(i) + residuals(i-1) * rhoo ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);     
      nll+=-dnorm(obs_logRS(i),pred_logRS(i),sigAR,true);      
    } 
  }
  
  //Type umsy     = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  //Type Smsy     = alpha/beta * (Type(0.5) -Type(0.07) * alpha);

  Type umsy = (Type(1) - LambertW(exp(1-alpha)));
  Type Smsy = (Type(1) - LambertW(exp(1-alpha))) / beta;
  
  Type ans = nll + pnll;

  REPORT(alpha)
  REPORT(beta)
  REPORT(rhoo)
  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(sigobs)
  REPORT(sigAR)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(nll);
  REPORT(pnll);  
 
  ADREPORT(alpha);
  ADREPORT(beta);
  REPORT(rhoo);
  ADREPORT(sigobs);
  ADREPORT(sigAR);
  ADREPORT(umsy);
  ADREPORT(Smsy);
  

  return ans;
}

