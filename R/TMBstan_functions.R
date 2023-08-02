# Functions to run the TMB stock recruitment  models
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Simple Ricker model estimated with TMBstan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. See details for priors documentation.
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param chains number of MCMC chains
#' @param iter number of MCMC iterations per chain
#' 
#'
#'
#' @details Priors: Weakly informative priors are included for the main parameterst of the model:
#' alpha ~ gamma(3,1)
#' logbeta ~ N(-12,3)
#' sigobs ~ gamma(2,1/3) 
#' 
#' 
#' @returns a list containing several model outputs - posterior distributions:
#' * alpha - posterior distribution for the alpha parameter vector
#' * beta - posterior distribution for the beta parameter 
#' * sigobs - posterior distribution for the observation error sigma   
#' * Smax - posterior distribution for the Smax parameter vector
#' * umsy - posterior distribution for the umsy parameter 
#' * Smsy - posterior distribution for the Smsy           
#' * AICc - AICc values, given by 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1)), excluding prior components
#' * BIC - BIC values, excluding prior components
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
ricker_TMBstan <- function(data,  silent = FALSE, control = TMBcontrol(), 
  tmb_map = list(), AC=FALSE, priors_flag=1, stan_flag=0,sig_p_sd=1,
  chains=6,iter=5000, warmup = floor(iter/2)) {

  
  tmb_data <- list(
    obs_S = data$S,
    obs_logRS = data$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    y_oos=mean(data$logRS),
    x_oos=mean(data$S)
  )
  
  magS <- log10_ceiling(max(data$S))
  initlm <- lm(logRS~S, data=data)
  tmb_random <- NULL

  if(!AC){
    tmb_params <- list(
      alpha   = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1)
    )
    lowlimit <- c(-4,-20,log(0.01))
    hightlimit <- c(20,-4,log(2))

    tmb_obj <- TMB::MakeADFun(data = tmb_data, 
                             parameters = tmb_params, 
                             map = tmb_map,
                             random = tmb_random, 
                             DLL = "Ricker_simple", 
                             silent = silent)
  
  }else{
    tmb_params <- list(
      alpha   = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1),
      rho=0
    )

    lowlimit <- c(-4,-20,log(0.01),-1)
    hightlimit <- c(20,-4,log(2),1)
   

    tmb_obj <- TMB::MakeADFun(data = tmb_data,
                              parameters = tmb_params, 
                              map = tmb_map,
                              random = tmb_random, 
                              DLL = "Ricker_autocorr", 
                              silent = silent)
  }




  tmb_mcmc <- tmbstan::tmbstan(tmb_obj, chains=chains,
              iter=iter, init="random",
              lower=lowlimit , upper=hightlimit,
               control = list(adapt_delta = 0.98),
               seed = 123,
               warmup = warmup)

    mc <- extract(tmb_mcmc, pars=names(tmb_obj$par),
              inc_warmup=FALSE, permuted=FALSE)
    

    fit_summary <- summary(tmb_mcmc)

    posterior <- as.matrix(tmb_mcmc)
  
    alpha <- NULL
    sigobs <- NULL  
    beta <- NULL
    Smax <- NULL
    umsy <- NULL 
    Smsy <- NULL  
    pred_logR<-matrix(NA,nrow=nrow(posterior),ncol=data$S)
    pred_logRS<-matrix(NA,nrow=nrow(posterior),ncol=data$S)



    for(i in 1:nrow(posterior)){
      r <- tmb_obj$report(posterior[i,-ncol(posterior)])

      alpha[i] <- r$alpha
      sigobs[i] <- r$sigobs  
      beta[i] <- r$beta
      Smax[i] <- r$Smax
      umsy[i] <- r$umsy  
      Smsy[i] <- r$Smsy  
      pred_logR[i,] <- r$pred_logR
      pred_logRS[i,] <- r$pred_logRS
      ll[i,] <-r$ll
      
    }

    #calculate WAIC
 
    elpd_1 <- apply(ll,2,log_mean_exp) #
    AICc  <- -2*sum(elpd_1) + 2*npar +( 2*npar*(npar+1)/(nrow(data)-npar-1))
    BIC  <- -2*sum(elpd_1) + npar*log(nrow(data))
  

  
  
  structure(list(
    alpha = alpha,
    beta = beta,
    Smax = Smax,
    sig = sigobs,  
    umsy = umsy,  
    Smsy = r$Smsy,  
    pred_logR = pred_logR,
    pred_logRS = pred_logRS,
    AICc = AICc,
    BIC=BIC,
    fit_summary =fit_summary,
    posterior=posterior))

}





#' Ricker model with random walk in a, b or both parameters with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param tv.par Which parameters should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control list of controls, as it would be passed to stan()
#' @param ini_param Optional. A list with initial parameter guesses. The list should contain: alphao (a number),
#' logbeta (a number), logsigobs (a number), logsiga (a number), and alpha ( a vector with the same length as the data). 
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. See details for priors documentation.
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param siga_p_sd sd for half normal prior on sigma for alpha random walk parameter. default is 1.
#' @param sigb_p_sd sd for half normal prior on sigma for beta random walk parameter. default is 1.
#' @param chains number of MCMC chains
#' @param iter number of MCMC iterations per chain
#' @param laplace Use the laplace approximation for random effects, default is FALSE. (TRUE not implemented)
#' 
#' @details Priors: Weakly informative priors are included for the main parameterst of the model:
#' alpha ~ gamma(3,1)
#' logbeta ~ N(-12,3)
#' sigobs ~ gamma(2,1/3)
#' siga ~ gamma(2,1/3) 
#' sigb ~ gamma(2,1/3)  
#' 
#' 
#' @returns a list containing several model outputs:
#' * alpha - MLE estimates for the alpha parameter vector
#' * beta - MLE estimates for the beta parameter 
#' * sig - MLE estimates for the observation error standard deviation     
#' * siga - MLE estimates for the process error (variation in alpha) standard deviation
#' * sigb - MLE estimates for the process error (variation in beta) standard deviation   
#' * AICc - AICc values, given by 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1)), excluding prior components
#' * BIC - BIC values, excluding prior components
#' * model - opt object, generated by the `stats::nlminb()`  
#' * tmb_data - data provided to TMB,
#' * tmb_params - parameter intial guesses,
#' * tmb_map - mapping indicating which parameters should be fixed or estimated,
#' * tmb_obj - obj object, generated by the `TMB::MakeADFun()` 
#' * gradients - final convergence gradient
#' * bad_eig - eigen value
#' * call - original function call
#' * sd_report - MLE estimates and sdt Error estimates for main parameter estimates
#' * class - name of cpp model
#' 
#' 
#' @export
#' @examples 
#' data(harck)
#' ricker_rwa_TMB(data=harck)
#' 
#' 
ricker_rw_TMBstan <- function(data, tv.par=c('a','b','both'), silent = FALSE, 
  control = list(adapt_delta = 0.98), ini_param=NULL, tmb_map = list(), priors_flag=1, stan_flag=1,
  sig_p_sd=1, siga_p_sd=1, sigb_p_sd=1, logb_p_mean=-12,logb_p_sd=3,
   chains=6,iter=10000 ,laplace=FALSE, warmup = floor(iter/2),...) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    obs_S = data$S,
    obs_logRS = data$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd, 
    logb_p_mean=logb_p_mean,
    logb_p_sd=logb_p_sd
  )

  if(is.null(ini_param)){
    magS <- log10_ceiling(max(data$S))
    initlm<-lm(logRS~S, data=data)
  }

  if(tv.par=="a"){

    tmb_data$siga_p_sd=siga_p_sd

    if(is.null(ini_param)){
      tmb_params <- list(alphao   = initlm$coefficients[[1]],
                   logbeta = ifelse(initlm$coefficients[[2]]>0,
                                   log(1/magS),
                                   log(-initlm$coefficients[[2]])),
                   logsigobs = log(.5),
                   logsiga = log(.5),
                   alpha = rep(1,length(tmb_data$obs_S)))
    }else{
      tmb_params <-ini_param
    }
    tmb_random <- "alpha"
    tmb_obj <- TMB::MakeADFun(data = tmb_data, 
                              parameters = tmb_params, 
                              map = tmb_map,
                              random = tmb_random, 
                              DLL = "Ricker_tva", 
                              silent = silent)

    lowlimit <- c(0.01,-20,log(0.01),log(0.01))
    hightlimit <- c(20,-4,log(2),log(2))

   clss <- "Ricker_tva"
   npar <- 4
   npar_all <- 4+(length(data$S)-1)

  }else if(tv.par=="b"){

    tmb_data$sigb_p_sd=sigb_p_sd

    if(is.null(ini_param)){
     
      tmb_params <- list(logbetao = ifelse(initlm$coefficients[[2]]>0,
                                           log(1/magS),
                                           log(-initlm$coefficients[[2]])),
                        alpha   = max(initlm$coefficients[[1]],.5),                 
                        logsigobs = log(.6),
                        logsigb = log(.2),
                        logbeta=rep(ifelse(initlm$coefficients[[2]]>0,
                                           log(1/max(data$S)),
                                           log(-initlm$coefficients[[2]])),
                                    length(data$S)))
    }else{
      tmb_params <-ini_param
    }

    tmb_random <- "logbeta"

    tmb_obj <- TMB::MakeADFun(data = tmb_data, 
      parameters = tmb_params, 
      map = tmb_map,
      random = tmb_random, 
      DLL = "Ricker_tvlogb", 
      silent = silent)

    lowlimit <- c(-20,0.01,log(0.01),log(0.01))
    hightlimit <- c(-4,20,log(2),log(2))

    clss <- "Ricker_tvlogb"
    npar <- 4
    npar_all <- 4+(length(data$S)-1)

  }else if(tv.par=="both"){

    tmb_data$siga_p_sd=siga_p_sd
    tmb_data$sigb_p_sd=sigb_p_sd
    
    if(is.null(ini_param)){    
      tmb_params <- list(logbetao = ifelse(initlm$coefficients[[2]]>0,
                                           log(1/magS),
                                           log(-initlm$coefficients[[2]])),
                        alphao   = max(initlm$coefficients[[1]],.5),                 
                        logsigobs = log(.6),
                        logsiga = log(.2),
                        logsigb = log(.2),
                        logbeta=log(rep(ifelse(initlm$coefficients[[2]]>0,
                                               log(magS),
                                               -initlm$coefficients[[2]]),
                                        length(data$S))),
                        alpha = rep(1,length(tmb_data$obs_S))      
     )

    }else{
      tmb_params <-ini_param
    }
    tmb_random <- c("logbeta", "alpha")

    tmb_obj <- TMB::MakeADFun(data = tmb_data, 
      parameters = tmb_params, 
      map = tmb_map,
      random = tmb_random, 
      DLL = "Ricker_tva_tvb", 
      silent = silent)
  
    lowlimit <- c(-20,0.01,log(0.0),log(0.01),log(0.01))
    hightlimit <- c(-4,20,log(2),log(2),log(2))
    
    clss <- "Ricker_tva_tvb"
    npar <- 5
    npar_all <- 5+(length(data$S)-1)*2

  }else{
    stop(paste("tv.par",tv.par,"not recognized."))
  }
 

  #===================================
  # TMBstan fit
  #===================================
  #chains=6
  #iter=5000
  tmb_mcmc <- tmbstan::tmbstan(tmb_obj, chains=chains,
              iter=iter, init="random",
              lower=lowlimit , upper=hightlimit,
               control = control,
               laplace=laplace,
                warmup = warmup )

   
  fit_summary <- summary(tmb_mcmc)

  posterior <- as.matrix(tmb_mcmc)
  
  
  postlist<-as.list(as.data.frame(t(posterior[,-ncol(posterior)])))
   

  r <- lapply(postlist,tmb_obj$report)
     
      
  beta <-sapply(r,"[[","beta")
  alpha  <- sapply(r,"[[","alpha")
  sigobs <- sapply(r,"[[","sigobs")
  Smax <- sapply(r,"[[","Smax")
  umsy <-sapply(r,"[[","umsy")
  Smsy <- sapply(r,"[[","Smsy")
  ll <- sapply(r,"[[","Smsy")
  pred_logR <- sapply(r,"[[","pred_logR")
  pred_logRS <- sapply(r,"[[","pred_logRS")
    
      
  elpd_1 <- apply(ll,1,log_mean_exp) #
  AICc  <- -2*sum(elpd_1) + 2*npar +( 2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- -2*sum(elpd_1) + npar*log(nrow(data))
  
       
      
      
  
  #todo add alpha, beta and sigma parameter esitimates
  structure(list(
    alpha = alpha,
    beta = beta,
    Smax = Smax,
    sig = sigobs,  
    umsy = umsy,  
    Smsy = r$Smsy,  
    pred_logR = pred_logR,
    pred_logRS = pred_logRS,
    AICc = AICc,
    BIC=BIC,
    fit_summary =fit_summary,
    posterior=posterior))

}






# END
#***********************************************************************************
