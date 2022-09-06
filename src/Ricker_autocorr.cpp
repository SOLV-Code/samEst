#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
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
  Type ans= Type(0);
  ans -=dnorm(alpha,Type(0.0),Type(2.5),true);
  ans -=dnorm(logbeta,Type(-12.0),Type(3.0),true); 
  //ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
  ans -= dgamma(sigobs,Type(2.0),Type(1.0)/Type(3.0),true);

  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps) ;
 
  pred_logRS(0) = alpha - beta * obs_S(0) ;
  pred_logR(0) = pred_logRS(0) + log(obs_S(0));
  residuals(0) = obs_logRS(0) - pred_logRS(0);
    
  ans+= -dnorm(obs_logRS(0),pred_logRS(0),sigobs,true);

  for(int i=1;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      pred_logRS(i) = alpha - beta * obs_S(i) + residuals(i-1) * rhoo ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);     
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sigAR,true);      
    } 
  }
  
  Type umsy     = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  Type Smsy     = alpha/beta * (Type(0.5) -Type(0.07) * alpha);  

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
  return ans;
}

