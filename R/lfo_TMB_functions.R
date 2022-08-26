#TMB lfo function
# Dan greenberg and modified by Catarina for the package




#' Simple Ricker model estimated with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' 
#' @param tv.par string. Which parameters are time-varying? Options are c('static','alpha','beta','both')
#' @param L starting point for LFO-CV (min. 10)
#' 
#' 
#' @returns vector of lfo by year
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
tmb_mod_lfo_cv=function(data,tv.par=c('static','alpha','beta','both'),L=10){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)

  if(L<10){
  	warning("L values < 10 are not recommended, results may be unstable")
  }
  if(sum(c("S","logRS") %in% names(data))!=2){
  	stop("data must contain 'S' and 'logRS', please revise names(data)") 
  }

  if(tv.par=='static'){
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tmb<- rickerTMB(data=df_past)
      rs_pred_1b=fit_past_tmb$alpha-fit_past_tmb$beta*df_oos$S[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tmb$sig))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(exact_elpds_1b)
  }else if(tv.par=='alpha'){
  
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_a_tmb<- ricker_rwa_TMB(data=df_past)
    
      rs_pred_1b=fit_past_tv_a_tmb$alpha[i]-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      rs_pred_3b=mean(fit_past_tv_a_tmb$alpha[(i-2):i])-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      rs_pred_5b=mean(fit_past_tv_a_tmb$alpha[(i-4):i])-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=exp(fit_past_tv_a_tmb$sig)))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=exp(fit_past_tv_a_tmb$sig)))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=exp(fit_past_tv_a_tmb$sig)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }else if(tv.par=='beta'){
  	#stop("not defined")
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_b_tmb<- ricker_rwb_TMB(df=df_past)
      
      rs_pred_1b=fit_past_tv_b_tmb$alpha-fit_past_tv_b_tmb$beta[i]*df_oos$S[i + 1]
      rs_pred_3b=fit_past_tv_b_tmb$alpha-mean(fit_past_tv_b_tmb$beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5b=fit_past_tv_b_tmb$alpha-mean(fit_past_tv_b_tmb$beta[(i-4):i])*df_oos$S[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=exp(fit_past_tv_b_tmb$sig)))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=exp(fit_past_tv_b_tmb$sig)))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=exp(fit_past_tv_b_tmb$sig)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }else if(tv.par=='both'){
  	
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_ab_tmb<- ricker_rwab_TMB(data=df_past)
      
      rs_pred_1b=fit_past_tv_ab_tmb$alpha[i]-fit_past_tv_ab_tmb$beta[i]*df_oos$S[i + 1]
      rs_pred_3b=mean(fit_past_tv_ab_tmb$alpha[(i-2):i])-mean(fit_past_tv_ab_tmb$beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5b=mean(fit_past_tv_ab_tmb$alpha[(i-4):i])-mean(fit_past_tv_ab_tmb$beta[(i-4):i])*df_oos$S[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_1b,sd=exp(fit_past_tv_ab_tmb$sig)))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_3b,sd=exp(fit_past_tv_ab_tmb$sig)))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_5b,sd=exp(fit_past_tv_ab_tmb$sig)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }else if(tv.par=='HMM'){
    #stop("not defined")
  }else{
  	stop(paste("tv.par", tv.par,"not defined, valid options are c('static','alpha','beta','both')"))
  }
}

