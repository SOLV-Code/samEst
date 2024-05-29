# Functions to run the RTMB stock recruitment  models
#===============================================









#' Simple Ricker model function to be called with RTMB
#'
#' @param par A list or data frame containing parameter names and initial guesses 
#' 
#'
#' @details  this function will
#
#' 
#' @returns penalized (if priors_flag= TRUE) negative log likelihood
#' 
#' 
#' @export
#' 
#' 
#' 
ricker_RTMB_fn <- function(param,dat){

  RTMB::getAll(param, dat)
  
  #parameters
  beta <- exp(logbeta)
  sigobs <- exp(logsigobs);
  Smax  <- 1/beta;

  if(priors_flag == 1){
    
    pnll <- -dnorm(logalpha, mean=1.5, sd=2.5, log=TRUE)
    pnll <- pnll - dnorm(logbeta,mean=logb_p_mean, sd=logb_p_sd, log=TRUE)
    pnll <- pnll - dnorm(sigobs,mean=0, sd=sig_p_sd, log=TRUE) -  pnorm(0, 0,sd=sig_p_sd, log.p=TRUE)

  }
  if(stan_flag){
    pnll <- pnll - logsigobs #Jacobian for half normal prior
  }
  
  pred_logRS <- logalpha - beta * obs_S 
  nll <- - sum(dnorm(obs_logRS, mean=pred_logRS, sd=sigobs, log=TRUE))

  pred_logR <- pred_logRS + log(obs_S)
  residuals <- obs_logRS - pred_logRS
    
  #RTMB does not accept the lambertW function
  Smsy <- (1 - LambertW0(exp(1 - logalpha))) /beta
  umsy <- (1 - LambertW0(exp(1 - logalpha)))
  Sgen <-  -1/beta*LambertW0(-beta*Smsy/exp(logalpha))
  
  
  ans <- nll + pnll;

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(logalpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Sgen)
  REPORT(residuals)
  REPORT(nll)
  REPORT(pnll)



  return(ans)
   

}









#' Simple Ricker model function to be called with RTMB
#'
#' @param par A list or data frame containing parameter names and initial guesses 
#' 
#'
#' @details  this function will
#
#' 
#' @returns penalized (if priors_flag= TRUE) negative log likelihood
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' 
#' 
ricker_ac_RTMB_fn <- function(param,dat){

  RTMB::getAll(param, dat)

  #Indexing
  N <- length(by)

  #transformed parameters
  beta <- exp(logbeta)
  sigobs <- exp(logsigobs);
  Smax  <- 1/beta;
  rho <- minus_one_to_one(ar1_phi)
  sigAR  <- sigobs*sqrt(1-rho^2)
  
  #empty objects
  residuals<-numeric(N)

  
  if(priors_flag == 1){
    
    pnll <- -dnorm(logalpha, mean=1.5, sd=2.5, log=TRUE)
    pnll <- pnll - dnorm(logbeta,mean=logb_p_mean, sd=logb_p_sd, log=TRUE)
    pnll <- pnll - dnorm(sigobs,mean=0, sd=sig_p_sd, log=TRUE) -  pnorm(0, 0,sd=sig_p_sd, log.p=TRUE)
    pnll <- pnll - dnorm(rho,mean=0, sd=1, log=TRUE) 
  }
  if(stan_flag){
    pnll <- pnll - logsigobs #Jacobian for half normal prior
    pnll <- pnll - (log(2.) + ar1_phi - 2. * log(1. + exp(ar1_phi))) #Jacobian for minus_one_to_one
  }
  
  pred_logRS <- logalpha - beta * obs_S 
  residuals[1] <- obs_logRS[1] - pred_logRS[1]
  nll <- -dnorm(obs_logRS[1], mean=pred_logRS[1], sd=sigobs, log=TRUE)
  
  for(i in 2:N){
    pred_logRS[i] <- pred_logRS[i] + (rho^(by[i]-by[i-1])*residuals[i-1])
    nll <- nll - dnorm(obs_logRS[i], mean=pred_logRS[i], sd=sigAR, log=TRUE)
    residuals[i] <- obs_logRS[i] - pred_logRS[i]
  }

  pred_logR <- pred_logRS + log(obs_S)
  
    
  
  #RTMB does not accept the lambertW function
  Smsy <- (1 - LambertW0(exp(1 - logalpha))) /beta
  umsy <- (1 - LambertW0(exp(1 - logalpha)))
  Sgen <-  -1/beta*LambertW0(-beta*Smsy/exp(logalpha))
  
  ans <- nll + pnll;

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(logalpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(rho)
  REPORT(sigAR)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Sgen)
  REPORT(residuals)
  REPORT(nll)
  REPORT(pnll)

  return(ans)
   

}















