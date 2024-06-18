#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_IVECTOR( options_z );
  DATA_VECTOR( obs_logRS );
  DATA_UPDATE( obs_logRS );
  DATA_VECTOR( obs_S );

  // Parameters
  PARAMETER( alpha );
  PARAMETER( beta );
  PARAMETER( ln_sigA );
  PARAMETER( ln_sigB );
  PARAMETER( ln_sigma );
  PARAMETER_VECTOR( epsA_t );
  PARAMETER_VECTOR( epsB_t );

  Type neglogpen = 0;
  // Define distribution for random effects
  for( int t=0; t<obs_logRS.size(); t++ ){
    if( options_z(0)==1 ){
      if(t==0) neglogpen -= dnorm( epsA_t(t), Type(0.0), exp(ln_sigA), true );
      if(t>=1) neglogpen -= dnorm( epsA_t(t), epsA_t(t-1), exp(ln_sigA), true );
    }
    if( options_z(1)==1 ){
      if(t==0) neglogpen -= dnorm( epsB_t(t), Type(0.0), exp(ln_sigB), true );
      if(t>=1) neglogpen -= dnorm( epsB_t(t), epsB_t(t-1), exp(ln_sigB), true );
    }
  }

  // Define distribution for data
  vector<Type> yhat_t( obs_logRS.size() );
  vector<Type> loglik_t( obs_logRS.size() );
  for( int t=0; t<obs_logRS.size(); t++ ){
    yhat_t(t) = (alpha + epsA_t(t)) + (beta + epsB_t(t)) * obs_S(t);
    loglik_t(t) = dnorm( obs_logRS(t), yhat_t(t), exp(ln_sigma), true );
  }

  // Report predictions
  Type jnll = neglogpen - sum(loglik_t);
  REPORT( yhat_t );
  REPORT( jnll );
  REPORT( loglik_t );
  ADREPORT( yhat_t );
  return jnll;
}
