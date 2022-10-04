#include <TMB.hpp>


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

template<class Type>
Type objective_function<Type>::operator() ()
{
  //Kalman filter translated to TMB from Carrie holt R-code
  DATA_VECTOR(x);   // observed log recruitment
  DATA_VECTOR(y);    // observed  Spawner
  
 
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //logsige     -> Natural log of the standard deviation of observation error
  //logsigw     -> Natural log of the standard deviation of system error

  PARAMETER(initmeana);
  //PARAMETER(loginitvara);
  PARAMETER(b);
  PARAMETER(logsige);
  PARAMETER(logsigw);

  Type sige = exp(logsige);
  Type sigw = exp(logsigw);
  Type sig;
  //Type initvara = exp(loginitvara);
  Type initvara = Type(1.0);

  int Tmax = x.size();

  vector<Type> priormeana(Tmax),  priorvara(Tmax), yhat(Tmax),
  f(Tmax), v(Tmax), postmeana(Tmax), postvara(Tmax), filtery(Tmax),
   smoothemeana(Tmax), smoothevara(Tmax), smoothey(Tmax), dummy(Tmax);

  vector<Type> umsy(Tmax), Smsy(Tmax);

  vector<Type> pstar(Tmax-1);
  
  Type ans = Type(0.0);

  
  
  //first year
  priormeana(0) = initmeana;
  priorvara(0) = initvara;
  yhat(0) = priormeana(0) + b * x(0);
  v(0) = y(0) - yhat(0);
  f(0) = priorvara(0) + sige*sige;
     
  postmeana(0) = priormeana(0) + (priorvara(0) * (v(0)/f(0)));
  postvara(0) = priorvara(0) - (priorvara(0)*priorvara(0)/f(0));
  filtery(0) = postmeana(0) + b * x(0);
  sig = sqrt(f(0));
  ans -=  dnorm(y(0) , yhat(0),sig,true);

  //ans +=  (log(f(0)) + (v(0)*v(0)/f(0)))/Type(2.);



  for(int t=1; t<Tmax; t++){

    
    priormeana(t) = postmeana(t-1);
    priorvara(t) = postvara(t-1) + sigw*sigw;

    if(isNA(x(t))){
      postmeana(t) = priormeana(t);
      postvara(t) = priorvara(t);
    }
     
    
    //Step 2: Generate predicted value for y(t) given y(t-1) and error
    if(!isNA(x(t))){
      yhat(t) = priormeana(t) + b * x(t);
      v(t) = y(t) - yhat(t);
      f(t) = priorvara(t) + sige*sige;
      
      // Step 3: Generate posterior distribution for intercept (a):
      postmeana(t) = priormeana(t) + (priorvara(t) * (v(t)/f(t)));
      postvara(t) = priorvara(t) - (priorvara(t)*priorvara(t)/f(t));
      filtery(t) = postmeana(t) + b * x(t);
      sig = sqrt(f(t));
      ans -=  dnorm(y(t) , yhat(t),sig,true);
      //ans += (log(f(t)) + (v(t)*v(t)/f(t)))/Type(2.);
    }

  }
  //End loop over time


  //smoothing part
  // Step 6: Smoothing of kalman filter estimates for time-varying intercept
  // Start loop over time (NB: Calculations start with last values first)
  
  smoothemeana(Tmax-1) = postmeana(Tmax-1);
  smoothevara(Tmax-1) = postvara(Tmax-1);
  smoothey(Tmax-1) = smoothemeana(Tmax-1) + b * x(Tmax-1);
  pstar.setZero();

  for(int i=Tmax-2; i>=0; --i){
      
      pstar(i) = postvara(i)/priorvara(i + 1);
      smoothemeana(i) = postmeana(i) + pstar(i) * (smoothemeana(i + 1) - priormeana(i + 1));
      smoothevara(i) = postvara(i) + pstar(i)*pstar(i) * (smoothevara(i + 1) - priorvara(i + 1));    
      smoothey(i) = smoothemeana(i) + b * x(i);
  }

  
  for(int i=0; i<Tmax; i++){
    Smsy(i) = (Type(1.0) - LambertW(exp(1-smoothemeana(i))) ) / -b;
    umsy(i) = (Type(1.0) - LambertW(exp(1-smoothemeana(i))) ); 
  }
  
   REPORT(ans)
   REPORT(Tmax)
   REPORT(priormeana)
   REPORT(priorvara)
   REPORT(yhat)
   REPORT(f)
   REPORT(v)
   REPORT(postmeana)
   REPORT(postvara)
   REPORT(filtery)
   REPORT(pstar)
   REPORT(smoothemeana)
   REPORT(smoothevara)
   REPORT(sige)
   REPORT(sigw)
   REPORT(smoothey)
   REPORT(initmeana)
   REPORT(initvara)
   REPORT(b)
   REPORT(Smsy)
   REPORT(umsy)

   ADREPORT(smoothemeana)
   ADREPORT(smoothevara)
   ADREPORT(sige)
   ADREPORT(sigw)
    return ans;

}

