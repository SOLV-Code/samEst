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
  DATA_SCALAR(sigb_p_sd); //sd for sigb prior

  DATA_SCALAR(logb_p_mean); //mean for logb prior
  DATA_SCALAR(logb_p_sd); //sd for logb prior

  //DATA_SCALAR(prbeta1); //beta prior parameter
  //DATA_SCALAR(prbeta2); //beta prior parameter
  
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(logbetao);
  PARAMETER(alpha);
  //PARAMETER(rho);
  //PARAMETER(logvarphi);
  PARAMETER(logsigobs);
  PARAMETER(logsigb);


  PARAMETER_VECTOR(logbeta);
  
  
  int timeSteps=obs_logRS.size();
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error (beta) sd
  
  //Type varphi     = exp(logvarphi);
  //Type theta     = sqrt(Type(1.0)/varphi);
  //Type sig       = sqrt(rho) * theta;
  //Type tau        = sqrt(Type(1.0)-rho) * theta ;

  Type sigobs = exp(logsigobs);
  Type sigb = exp(logsigb);
  

  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps), Smsy(timeSteps), residuals(timeSteps);
  vector<Type> Srep(timeSteps), beta(timeSteps), Smax(timeSteps), ll(timeSteps);

  

  //prior on observation and process variance ratio
  //Type ans= -dbeta(rho,prbeta1,prbeta2,true);  
  //priors on parameters
  //Type ans= Type(0);
  Type nll = Type(0.0);
  Type pnll = Type(0.0);
  Type renll = Type(0.0);

  if(priors_flag == 1){
    //ans -= dnorm(alpha,Type(0.0),Type(2.5),true);
    pnll -=dnorm(alpha,Type(1.5),Type(2.5),true);
    //pnll -= dgamma(alpha,Type(3.0),Type(1.5),true);
    //pnll -=dnorm(alpha,Type(2.5),Type(3.0),true)- log(pnorm(Type(-2.0), Type(2.5),Type(3.0)));
    pnll -= dnorm(logbetao,logb_p_mean,logb_p_sd,true);
    //ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
    //ans -= dnorm(logsigb,Type(0.0),Type(2.0),true);
    //ans -= dnorm(sigobs,Type(0.0),Type(2.0),true);
    //ans -= dnorm(sigb,Type(0.0),Type(2.0),true);
    //pnll -= dgamma(sigobs,Type(2.0),Type(1.0)/Type(3.0),true);
    //pnll -= dgamma(sigb,Type(2.0),Type(1.0)/Type(3.0),true);
    
    pnll  -= dnorm(sigobs,Type(0.0),sig_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sig_p_sd));
    pnll  -= dnorm(sigb,Type(0.0),sigb_p_sd,true) - log(pnorm(Type(0.0), Type(0.0),sigb_p_sd));
   
    if(stan_flag){
      pnll -= logsigobs;
      pnll -= logsigb;
    } 

  }
  
  renll+= -dnorm(logbeta(0),logbetao,sigb,true);
  
  
  for(int i=1;i<timeSteps;i++){  
    renll+= -dnorm(logbeta(i),logbeta(i-1),sigb,true);  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      beta(i) = exp(logbeta(i));
      Smax(i) = Type(1.0)/beta(i);
      pred_logRS(i) = alpha - beta(i) * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      
      // Use the Hilborn approximations for Smsy and umsy
      Srep(i) = alpha/beta(i);
      //Smsy(i) =  alpha/beta(i) * (Type(0.5) -Type(0.07) * alpha);
      Smsy(i) = (1 - LambertW(exp(Type(1)-alpha)) ) / beta(i);

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ll(i)= dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
      nll+=-ll(i);
    }
  }
   // Use the Hilborn approximations for Smsy and umsy
  //Type umsy  = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  Type umsy  = (Type(1) - LambertW(exp(1-alpha)) ); 
 
  Type ans = nll + renll + pnll;

  REPORT(pred_logRS)
  REPORT(alpha)
  REPORT(sigobs)
  REPORT(sigb)
  REPORT(residuals)
  REPORT(beta)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)
  REPORT(nll);
  REPORT(ll);
  REPORT(pnll); 
  REPORT(renll); 

  ADREPORT(alpha);
  ADREPORT(beta);
  ADREPORT(sigobs);
  ADREPORT(umsy);
  ADREPORT(Smsy);
  return ans;
}

