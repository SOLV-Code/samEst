
#=====================================================
#RTMB

#data<-data.frame(by=c(1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016),
#    S=c(119511,173992,162077, 78760, 34800, 74145,176780, 90254,129841,118390, 98044, 28490, 39712, 71304,202362,105219, 84607,114067, 90448,247401,141958, 88801, 60786, 81122, 43156, 76459,106317,125251, 47664, 47935, 47075,102334, 46158),
#    R=c(211645, 82480,620514,245022,398430,311483,479719, 75985, 59537, 68779,326082, 80104,210970,147508,277290,497697,276425,154652, 86189,196445, 58069,133144, 63883,241104, 70693, 72491, 70075,141107, 45729, 44591, 85305, 72173, 67202),
#    logRS=c(0.57150193,-0.74645348, 1.34247664, 1.13494275, 2.43791444, 1.43532213, 0.99829449,-0.17209197,-0.77971266,-0.54308579, 1.20173253, 1.03377264, 1.67006253, 0.72692998, 0.31500572, 1.55394757, 1.18392253, 0.30439144,-0.04823254,-0.23062801,-0.89389928, 0.40503334, 0.04969379, 1.08927418, 0.49352510,-0.05329223,-0.41685910, 0.11919874,-0.04144373,-0.07231388, 0.59449099,-0.34917596, 0.37563272))


library(RTMB)
devtools::document()
devtools::load_all()

#============================
#simple ricker
p <- ricker_TMB(data=harck,priors_flag=0)
p2 <- ricker_RTMB(data=harck,priors_flag=0)

p$alpha
p2$logalpha

p$Smax
p2$Smax

p$tmb_obj$fn()
p2$obj$fn()



pac <- ricker_TMB(data=harck,AC=TRUE,priors_flag=0)
#This does not work
pac2 <- ricker_RTMB(data=harck,AC=TRUE,priors_flag=0)
#Error in dnorm(residuals[1], 0, sigobs, log = TRUE) : 
#  Non-numeric argument to mathematical function
#this works but produces different estimates than pac
pac3 <- ricker_RTMB(data=harck,AC=TRUE,,priors_flag=0,dautoreg=1)

#TMB with RTMB style of residuals in nll
pac4 <- ricker_TMB(data=harck,AC=TRUE,priors_flag=0,rtmb_v=1)


pac$rho
#these are the same but I cannot tell why the results are so different. 
pac2$rho
pac3$rho
pac4$rho



pac$alpha
#pac2$logalpha
pac3$logalpha

#likelihoods are similar though
pac$tmb_obj$fn()
pac2$obj$fn()
pac3$obj$fn()



#run AR models locally - outside of package
ricker_ac_RTMB_fn_local <- function(param,dat){

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
  pred_logRS <-numeric(N)
  residuals <-numeric(N)
  pred_logR <-numeric(N)
  pnll <- 0
  nll <- 0 
  
  if(priors_flag == 1){
    
    pnll <- -dnorm(logalpha, mean=1.5, sd=2.5, log=TRUE)
    pnll <- pnll - dnorm(logbeta, mean=logb_p_mean, sd=logb_p_sd, log=TRUE)
    pnll <- pnll - dnorm(sigobs, mean=0, sd=sig_p_sd, log=TRUE) -  pnorm(0, 0,sd=sig_p_sd, log.p=TRUE)
    pnll <- pnll - dnorm(rho, mean=0, sd=1, log=TRUE) 

    if(stan_flag == 1){
      pnll <- pnll - logsigobs #Jacobian for half normal prior
      pnll <- pnll - (log(2.) + ar1_phi - 2. * log(1. + exp(ar1_phi))) #Jacobian for minus_one_to_one
    }
  }
  
  pred_logRS <- logalpha - beta * obs_S

  residuals <- obs_logRS - pred_logRS
  
  
  nll <- nll - dnorm( residuals[1],0.0,sigobs,log=TRUE)
  for(i in 2:N){
    nll <- nll - dnorm(residuals[i],rho*residuals[i-1],sigAR,log=TRUE)
  }


  ans <- nll + pnll;

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(logalpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(rho)
  REPORT(sigAR)
  REPORT(sigobs)
  REPORT(Smax)
 
  REPORT(residuals)
  REPORT(nll)
  REPORT(pnll)

  return(ans)
   

}



data=harck
dat<- list(
    obs_S = data$S,
    obs_logRS = data$logRS,
    by = data$by,
    priors_flag=0,
    stan_flag=0,
    sig_p_sd=1,
    logb_p_mean=-12,
    logb_p_sd=3
  )
  

  magS <- log10_ceiling(max(data$S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

 
param <- list(
    logalpha   = initlm$coefficients[[1]],
    logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
    logsigobs = log(1),
    ar1_phi=0)



ricker_ac_RTMB_fn_local(param,dat)
ar1tmbcpfn<-function(param){ricker_ac_RTMB_fn_local(param,dat=dat)}
obj <- MakeADFun(ar1tmbfn, param, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
#similar to adautoreg results
obj$report()$rho


#==== 
#direct translation of TMB model

ricker_ac_RTMB_fn_local <- function(param,dat){

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
  pred_logRS <-numeric(N)
  residuals <-numeric(N)
  pred_logR <-numeric(N)
  pnll <- 0
  nll <- 0 
  
  if(priors_flag == 1){
    
    pnll <- -dnorm(logalpha, mean=1.5, sd=2.5, log=TRUE)
    pnll <- pnll - dnorm(logbeta, mean=logb_p_mean, sd=logb_p_sd, log=TRUE)
    pnll <- pnll - dnorm(sigobs, mean=0, sd=sig_p_sd, log=TRUE) -  pnorm(0, 0,sd=sig_p_sd, log.p=TRUE)
    pnll <- pnll - dnorm(rho, mean=0, sd=1, log=TRUE) 

    if(stan_flag == 1){
      pnll <- pnll - logsigobs #Jacobian for half normal prior
      pnll <- pnll - (log(2.) + ar1_phi - 2. * log(1. + exp(ar1_phi))) #Jacobian for minus_one_to_one
    }
  }
  
  pred_logRS[1] <- logalpha - beta * as.numeric(obs_S[1])
  residuals[1] <- obs_logRS[1] - pred_logRS[1]
  nll <- nll - dnorm(obs_logRS[1], pred_logRS[1], sigobs,log=TRUE)
  
  print(class(obs_logRS[1]))
  print(class(pred_logRS[1]))
  print(class(sigobs))
  print(class(nll))

 for(i in 2:N){
    pred_logRS[i] <- logalpha - beta * obs_S[i] + residuals[i-1] * rho 
    residuals[i] <- obs_logRS[i] - pred_logRS[i]
    nll <- nll - dnorm(obs_logRS[i],pred_logRS[i],sigAR,log=TRUE)
  }


  ans <- nll + pnll;

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(logalpha)
  REPORT(beta)
  REPORT(sigobs)
  REPORT(rho)
  REPORT(sigAR)
  REPORT(sigobs)
  REPORT(Smax)
 
  REPORT(residuals)
  REPORT(nll)
  REPORT(pnll)

  return(ans)
   

}





ricker_ac_RTMB_fn_local(param,dat)
ar1tmbcopyfn<-function(param){ricker_ac_RTMB_fn_local(param,dat=dat)}
obj_tmbcp <- MakeADFun(ar1tmbcopyfn, param, silent=TRUE)
opt_tmbcp <- nlminb(obj_tmbcp$par, obj_tmbcp$fn, obj_tmbcp$gr)
#similar to adautoreg results
obj_tmbcp$report()$rho


