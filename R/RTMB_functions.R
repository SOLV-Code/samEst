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
ricker_ac_RTMB_fn <- function(par,dat){

  RTMB::getAll(param, dat)

  #parameters
  beta <- exp(logbeta)
  sigobs <- exp(logsigobs);
  Smax  <- 1/beta;
  rhoo <- minus_one_to_one(rho)
  sigAR  <- sigobs*sqrt(1-rhoo^2)

  #dims
  N <-length(obs_S) #lnumber of observations
  ii <- rep(1,L)
  ii[is.na(obs_logRS)]<-0


  if(priors_flag == 1){
    
    pnll <- -dnorm(logalpha, mean=1.5, sd=2.5, log=TRUE)
    pnll <- pnll - dnorm(logbeta,mean=logb_p_mean, sd=logb_p_sd, log=TRUE)
    pnll <- pnll - dnorm(sigobs,mean=0, sd=sig_p_sd, log=TRUE) -  pnorm(0, 0,sd=sig_p_sd, log.p=TRUE)

  }
  if(stan_flag){
    pnll <- pnll - logsigobs #Jacobian for half normal prior
  }
  
  pred_logRS <- logalpha - beta * obs_S 
  residuals <- obs_logRS- pred_logRS


  nll <- -dnorm(obs_logRS[1], mean=pred_logRS[1], sd=sigobs, log=TRUE)
  for(i in 2:N){
    #pred_logRS[i] <- pred_logRS[i] + (rho^(ii[i]-ii[i-1])*epsilon[i-1])
    #nll <- nll - dnorm(obs_logRS[i], mean=pred_logRS[i], sd=sigobs, log=TRUE))
  }

  pred_logR <- pred_logRS + log(obs_S)
  residuals <- obs_logRS - pred_logRS
    
  nll <- -dnorm(obs_logRS[1], mean=pred_logRS[1], sd=sigobs, log=TRUE)

  #RTMB does not accept the lambertW function
  #Smsy <- (1 - lamW(exp(1 - logalpha))) /beta
  #umsy <- (1 - samest_lambertW(exp(1 - logalpha)))
  
  
  ans <- nll + pnll;

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(logalpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(Smax)
  #REPORT(umsy)
  #REPORT(Smsy)
  REPORT(residuals)
  REPORT(nll)
  REPORT(pnll)



  return(ans)
   

}








#' Simple Ricker model estimated with RTMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. 
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param logb_p_mean mean for prior on log b, default is -12.
#' @param logb_p_sd sd for prior on log b, default is 3.#' 
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
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
ricker_RTMB <- function(data,  silent = FALSE, control = TMBcontrol(), 
  tmb_map = list(), AC=FALSE, priors_flag=1, stan_flag=0,sig_p_sd=1,
  logb_p_mean=-12,logb_p_sd=3){

  
  dat<- list(
    obs_S = data$S,
    obs_logRS = data$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logb_p_mean=logb_p_mean,
    logb_p_sd=logb_p_sd
  )
  
  magS <- log10_ceiling(max(dat$obs_S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

  tmbfn<-function(param){ricker_RTMB_fn(param,dat=dat)}


  if(!AC){
    param <- list(
      logalpha  = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1)
    )

    obj <- RTMB::MakeADFun(tmbfn, parameters=param, 
      map = tmb_map,
      silent=silent)
  
  }else{
    param <- list(
      logalpha   = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1),
      rho=0
    )

    obj <- RTMB::MakeADFun(rickerac_RTMB_fn, param, 
      map = tmb_map,
      silent=TRUE)
  }

  opt <- stats::nlminb(obj$par, 
                       obj$fn, 
                       obj$gr,
                       control = control)

  sd_report <- RTMB::sdreport(obj)
  conv <- get_convergence_diagnostics(sd_report)

  nll <- obj$fn()[1]
  npar <- length(param)
 
  AICc  <- 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- 2*nll + npar*log(nrow(data))


  
  structure(list(
    logalpha   = obj$report()$logalpha,
    beta       = obj$report()$beta,
    Smax       = obj$report()$Smax,
    sig        = obj$report()$sigobs,
    sigar      = ifelse(AC,tmb_obj$report()$sigAR,NA),
    rho        = ifelse(AC,tmb_obj$report()$rhoo,NA),
    Smsy        = obj$report()$Smsy,
    umsy        = obj$report()$umsy,
    Sgen        = obj$report()$Sgen,
    AICc       = AICc,
    BIC        = BIC,
    residuals  = obj$report()$residuals,
    model      = opt,
    data   = data,
    params = param,
    map    = tmb_map,
    obj    = obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    conv_problem= conv$conv_problem,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = ifelse(AC,"Ricker_autocorr","Ricker_simple"))


}


  








