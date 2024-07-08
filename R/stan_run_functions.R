#' Simple Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series.
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param smax_priors Options for custom Smax (capacity) priors - 2 values for the mean and the scale (variance). Defaults to mean = half maximum observed spawners, scale = max. observed spawners (from data) 
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' 
#' @param ... Anything else that would be passed to rstan::sampling
#' 
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples
#' 
#' 
#' @importFrom rstan stan extract summary 
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerstan(data=harck)
#' 
ricker_stan <- function(data,  AC=FALSE, smax_priors=NULL, control = stancontrol(), mod=NULL, warmup=300, chains = 6, iter = 1000,...) {
 
   sm=samEst::sr_mod(type='static',ac=AC,par='n')

  if(AC==TRUE){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=seq_len(nrow(data)),
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)/2)
    
    if(is.null(smax_priors)==FALSE){
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=seq_len(nrow(data)),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
    }
    
  }else{
    datm = list(N=nrow(data),
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)/2)
    if(is.null(smax_priors)==FALSE){
      datm = list(N=nrow(data),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
    }
  }
  
  fit<-rstan::sampling(sm, data=datm,
                      control = control, warmup = warmup, 
                      chains = chains, iter = iter,verbose=FALSE)
  #fit <- rstan::stan(model_code = sm, 
  #                      data = datm,
  #                      control = control, warmup = warmup, chains = chains, iter = iter)
  #
  if(AC==FALSE){
    mc <- rstan::extract(fit,pars=c('log_a','S_max','U_msy','S_msy','sigma'),
                         inc_warmup=FALSE, permuted=FALSE)
  }
  if(AC==TRUE){
    mc <- rstan::extract(fit,pars=c('log_a','S_max','U_msy','S_msy','sigma','rho'),
                         inc_warmup=FALSE, permuted=FALSE)
  }
    
  aa <- rstan::summary(fit)
  
  ans<-list(alpha=c(median=aa$summary["log_a","50%"],sd=aa$summary["log_a","sd"],est2.5=aa$summary["log_a","2.5%"],est97.5=aa$summary["log_a","97.5%"]),
   beta=c(median=aa$summary["b","50%"],sd=aa$summary["b","sd"],est2.5=aa$summary["b","2.5%"],est97.5=aa$summary["b","97.5%"]),
   Smax=c(median=aa$summary["S_max","50%"],sd=aa$summary["S_max","sd"],est2.5=aa$summary["S_max","2.5%"],est97.5=aa$summary["S_max","97.5%"]),
   Smsy=c(median=aa$summary["S_msy","50%"],sd=aa$summary["S_msy","sd"],est2.5=aa$summary["S_msy","2.5%"],est97.5=aa$summary["S_msy","97.5%"]),
   Umsy=c(median=aa$summary["U_msy","50%"],sd=aa$summary["U_msy","sd"],est2.5=aa$summary["U_msy","2.5%"],est97.5=aa$summary["U_msy","97.5%"]),
   sigma=c(median=aa$summary["sigma","50%"],sd=aa$summary["sigma","sd"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]),
   summary=fit,
   samples=mc) 
  
  return(ans)

}








#' random walks Ricker model estimated with stan
#'
#' @param data A data frame containing Spawners (S) and log(Recruits/Spawners) (R_S) time series. Use sr_format for the correct column names. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' @param sm_ext Default is null, external stan rw model. Implemented to speed up simulation evaluation, use with caution.
#' @param lambertW Logical, indicating if lambertW functions should be used in stan code, requires git installation from stan
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples
#' 
#' 
#' @importFrom rstan stan extract summary 
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_rw_stan(data=harck)
#' 
ricker_rw_stan <- function(data, par=c('a','b','both'),smax_priors=NULL,  control = stancontrol(), mod=NULL,
  warmup=300,  chains = 6, iter = 1000, lambertW=FALSE,...) {
  #par='b'

  if(is.null(mod)){
   sm <- compile_code(type='rw',ac=FALSE,par=par,lambertW=lambertW)
  }else{
   sm <- mod
  }
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='rw',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <-sm_ext
  #}

  
  fit <- rstan::sampling(sm, data = list(N=nrow(data),
                                    L=max(data$by)-min(data$by)+1,
                                    ii=seq_along(data$logRS),
                                    R_S =data$logRS,
                                    S=data$S,
                                    pSmax_mean=max(data$S)/2,
                                    pSmax_sig=max(data$S)),
                        control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
  if(is.null(smax_priors)==FALSE){
    
    fit <- rstan::sampling(sm, data = list(N=nrow(data),
                                           L=max(data$by)-min(data$by)+1,
                                           ii=seq_along(data$logRS),
                                           R_S =data$logRS,
                                           S=data$S,
                                           pSmax_mean=smax_priors[1],
                                           pSmax_sig=smax_priors[2]),
                           control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
    
    
  }
  mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    
  aa <- rstan::summary(fit)
  
  ans <-list(alpha=aa$summary[grep("log_a",row.names(aa$summary)),"50%"],
   beta=aa$summary[grep("^b",row.names(aa$summary)),"50%"],
   Smax=aa$summary[grep("S_max",row.names(aa$summary)),"50%"],
   sigobs=aa$summary["sigma","50%"],
   siga=ifelse(par=="a"|par=="both",aa$summary["sigma_a","50%"],NA),
   sigb=ifelse(par=="b"|par=="both",aa$summary["sigma_b","50%"],NA),
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mc )

  if(lambertW){
    ans$Smsy=aa$summary["S_msy","50%"]
    ans$umsy=aa$summary["U_msy","50%"]
  }

  return(ans)

}



#' random walks Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param k_regime number of regimes in the model, default is 2
#' @param smax_priors Options for custom Smax (capacity) priors - 2 values for the mean and the scale (variance). Defaults to mean = half maximum observed spawners, scale = max. observed spawners (from data) 
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param sm_ext Default is null, external stan hmm model. Implemented to speed up simulation evaluation, use with caution.
#' @param lambertW Logical, indicating if lambertW functions should be used in stan code, requires git installation from stan
#' @param dirichlet_prior k_regime x k_regime matrix. Prior for transition probability matrix, 
#' if NULL prior is set to matrix(1,nrow=k_regime,ncol=k_regime)
#' @param ... Anything else that would be passed to rstan::sampling
#'  
#'
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples 
#' 
#' 
#' @importFrom rstan stan extract summary 
#'
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_hmm_stan(data=harck)
#' 
ricker_hmm_stan <- function(data, par=c('a','b','both'), k_regime=2, smax_priors=NULL, dirichlet_prior=NULL,
  control = stancontrol(), warmup=300,  chains = 6, iter = 1000, mod=NULL,
  lambertW=FALSE,...) {
  #par='both'
  
  if(is.null(mod)){
   sm <- compile_code(type='hmm',ac=FALSE,par=par,lambertW=lambertW)
  }else{
   sm <-mod
  }
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='hmm',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <- sm_ext
  #}
  if(is.null(dirichlet_prior)){
    dirichlet_prior<-matrix(c(2,1,1,2),nrow=k_regime,ncol=k_regime)
  }else if(nrow(dirichlet_prior)!=k_regime |ncol(dirichlet_prior)!=k_regime){
    stop("dirichlet_prior should be a k_regime x k_regime matrix")
  }

  
  fit <- rstan::sampling(sm, 
                        data = list(N=nrow(data),
                                    R_S =data$logRS,
                                    S=data$S,
                                    K=k_regime,
                                    alpha_dirichlet=dirichlet_prior,
                                    pSmax_mean=max(data$S)/2,
                                    pSmax_sig=max(data$S)),
                        control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
  
  if(is.null(smax_priors)==FALSE){
    fit <- rstan::sampling(sm, 
                           data = list(N=nrow(data),
                                       R_S =data$logRS,
                                       S=data$S,
                                       K=k_regime,
                                       alpha_dirichlet=dirichlet_prior,
                                       pSmax_mean=smax_priors[1],
                                       pSmax_sig=smax_priors[2]),
                           control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
    
  }

  mc <- rstan::extract(fit, 
          inc_warmup=FALSE, permuted=FALSE)
    

  aa <- rstan::summary(fit)
  parts <-stan_regime_rps(m=fit,par=par)
  
  ans<-list(
   alpha_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_t,NA),
   alpha_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_wt,NA),
   Smax_regime=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$S_max_t,NA),
   Smax_wgt=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$S_max_wt,NA),
   alpha=ifelse(par=="a"|par=="both",list(aa$summary[grep("log_a\\[",row.names(aa$summary)),"50%"]),aa$summary["log_a","50%"])[[1]],
   beta=ifelse(par=="b"|par=="both",list(aa$summary[grep("b\\[",row.names(aa$summary)),"50%"]),aa$summary["b","50%"])[[1]],
   Smax=ifelse(par=="b"|par=="both",list(aa$summary[grep("S_max\\[",row.names(aa$summary)),"50%"]),aa$summary["S_max","50%"])[[1]],
   sigobs=aa$summary["sigma","50%"],
   pi=aa$summary[grep("pi1",row.names(aa$summary)),"50%"],
   A=aa$summary[grep("A",row.names(aa$summary)),"50%"],
   probregime =matrix(aa$summary[grep("gamma\\[",row.names(aa$summary)),"50%"],ncol=k_regime, byrow=T),
   regime = aa$summary[grep("^zstar",row.names(aa$summary)),"50%"],
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mc ) 

  #extract time-series of parameters


  if(lambertW){
    
    
    ans$Smsy_regime=parts$S_msy_t
    ans$Smsy_wgt=parts$S_msy_wt
    ans$umsy_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$U_msy_t,NA)
    ans$umsy_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$U_msy_wt,NA)
    ans$Smsy=aa$summary[grep("S_msy\\[",row.names(aa$summary)),"50%"]
    ans$umsy=ifelse(par=="a"|par=="both",list(aa$summary[grep("U_msy\\[",row.names(aa$summary)),"50%"]),aa$summary["U_msy","50%"])[[1]]
   
  }

   return(ans)

}







#' Sampling control options. 
#'
#' Any arguments to pass to [rstan::sampling].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [rstan::sampling].
#'
#' @export
stancontrol <- function(adapt_delta = 0.99,  ...) {
  list(adapt_delta = 0.99, ...)
}



