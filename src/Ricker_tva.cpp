#include <TMB.hpp>

// Code taken from Brooke
// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
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
  DATA_INTEGER(priors_flag); //flag indicating wether or not priors should be used
  DATA_INTEGER(stan_flag); //flag indicating wether or not tmbstan is used 
  DATA_SCALAR(sig_p_sd); //sd for sigma prior
  DATA_SCALAR(siga_p_sd); //sd for siga prior

  PARAMETER(alphao); //initial alpha value
  PARAMETER(logbeta); //log of beta from ricker curve
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
 
  //Type ans= Type(0);
  Type nll = Type(0.0);
  Type renll = Type(0.0);
  Type pnll = Type(0.0);
  Type pnlla = Type(0.0);
  Type pnllb = Type(0.0);
  Type pnllc = Type(0.0);
  Type pnlld = Type(0.0);

  if(priors_flag == 1){
    //ans -=dnorm(alphao,Type(0.0),Type(2.5),true);
    pnll -=dgamma(alphao,Type(3.0),Type(1.0),true);
    pnll -=dnorm(logbeta,Type(-12.0),Type(3.0),true); 
    //ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
    //ans -= dnorm(logsiga,Type(0.0),Type(2.0),true);
    //pnll  -= dgamma(sigobs,Type(2.0),Type(1.0)/Type(3.0),true);  
    //pnll  -= dgamma(siga,Type(2.0),Type(1.0)/Type(3.0),true);

    pnll  -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    pnll  -= dnorm(siga,Type(0.0),siga_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),siga_p_sd));
    
    if(stan_flag){
      pnll -= logsigobs;
      pnll -= logsiga;
    } 
  }
  
  renll+= -dnorm(alpha(0),alphao,siga,true);
  
  for(int i=1;i<timeSteps;i++){  
    renll+= -dnorm(alpha(i),alpha(i-1),siga,true);
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      
      pred_logRS(i) = alpha(i) - beta * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));

      //umsy(i) = Type(.5) * alpha(i) - Type(0.07) * (alpha(i) * alpha(i)); 
      //Smsy(i) =  alpha(i)/beta * (Type(0.5) -Type(0.07) * alpha(i));
      Srep(i) = alpha(i)/beta;

      Smsy(i) = (Type(1) - LambertW(exp(1-alpha(i))) ) / beta;
      umsy(i) = (Type(1) - LambertW(exp(1-alpha(i))) ); 

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      nll+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
    }
  }

  Type ans= nll + renll + pnll;

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
  REPORT(nll);
  REPORT(pnll);  
  REPORT(pnlla);  
  REPORT(pnllb);  
  REPORT(pnllc);  
  REPORT(pnlld);  

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sigobs);
  ADREPORT(umsy);
  ADREPORT(Smsy);

  return ans;
}

