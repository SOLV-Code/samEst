#include <TMB.hpp>

// Code taken from Brooke
// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function


template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}


template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(obs_S);    // observed  Spawner
  DATA_VECTOR(obs_logRS);   // observed recruitment
  
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value

  PARAMETER(alphao);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  PARAMETER(logsiga);
 
  PARAMETER_VECTOR(alpha);
    
  int timeSteps=obs_logRS.size();

  Type beta = exp(logbeta);
  Type Smax = Type(1.0)/beta;
 
  Type sigobs = exp(logsigobs);
  Type siga = exp(logsiga);


  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps),umsy(timeSteps), Smsy(timeSteps), residuals(timeSteps);
  vector<Type> Srep(timeSteps); //, alpha(timeSteps)

 
  //priors on parameters
 
  Type ans= Type(0);
  ans -=dnorm(alphao,Type(0.0),Type(2.5),true);
  ans -=dnorm(logbeta,Type(-12.0),Type(3.0),true);
  
  ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
  ans -= dnorm(logsiga,Type(0.0),Type(2.0),true);

  //ans -= dgamma(logsigobs,Type(2.0),Type(0.333),true);
  //ans -= dgamma(logsiga,Type(2.0),Type(0.333),true);
  
  ans+= -dnorm(alpha(0),alphao,siga,true);
  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),siga,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      
      pred_logRS(i) = alpha(i) - beta * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));

      umsy(i) = Type(.5) * alpha(i) - Type(0.07) * (alpha(i) * alpha(i)); 
      Smsy(i) =  alpha(i)/beta * (Type(0.5) -Type(0.07) * alpha(i));
      Srep(i) = alpha(i)/beta;

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
    }
  
  }

  REPORT(alpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(siga)
  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(alphao)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)

  return ans;
}

