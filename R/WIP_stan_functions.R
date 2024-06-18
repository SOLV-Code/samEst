#' stan compile - compiles stan versions of static, autocorrelated, random walk, and regime shift models
#'
#' @param type Options for:
#' * all (default) - compiles 8 model sets
#' * static - compiles the classic (static) Ricker model
#' * ac - compiles the static Ricker model with 1-y residual autocorrelation
#' * rw - random walk (time-varying) model, must also specify the parameter ('par') - see below
#' * hmm - hidden markov (regime shift) model, must also specify the parameter ('par') - see below
#' @param par Options for:
#' * prod - sets stock productivity to vary through time or regime
#' * cap - sets stock capacity (Smax) to vary through time or regime
#' * both - sets both stock productivity and capacity to vary through time or regime (generally not recommended) 
#' @returns Nothing. It will compile stan code so these models can be sampled, to be used before 'stan_fit' (see help)
#' 
#' 
#' @import cmdstanr
#' @examples
#' stan_compile()
#' stan_compile(type='rw',par='prod')

stan_compile<-function(type=c('all'),par=c('prod','cap','both')) {
  if(type=='all'){
   m1=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[1])
   m2=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[2])
   m3=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[3])
   m4=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[4])
   m5=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[5])
   m6=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[6])
   m7=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[7])
   m8=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[8])
  }else if(type=='static'){
    m1=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[1])
  }else if(type=='ac'){
    m2=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[2])
  }else if(type=='rw'&par=='prod'){
    m3=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[3])
  }else if(type=='rw'&par=='cap'){
    m4=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[4])
  }else if(type=='rw'&par=='both'){
    m5=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[5])
  }else if(type=='hmm'&par=='prod'){
    m6=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[6])
  }else if(type=='hmm'&par=='cap'){
    m7=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[7])
  }else if(type=='hmm'&par=='both'){
    m8=cmdstanr::cmdstan_model(file.path('./src/stan/',list.files('./src/stan/','.stan'))[8])
  }else{
    print('Warning: model choice not found. type options: all, static, ac, rw, hmm; par options: a, b, both')
  }
}

#' stan_format - formats data for stan models
#'
#' @param s Time-series of spawner abundance
#' @param r Time-series of recruitment
#' @param t Time-series of brood years
#' @param K number of regimes for the hidden markov model - dfeault is 2 (high and low productivity or capacity)
#' @param smax_priors Custom priors for Smax - these should be two values (e.g. smax_priors=c(1000, 1000)) corresponding to the prior mean and prior error term.
#' @returns A list with all necessary inputs for the 8 different stan models
#' 
#' @examples
#' stan_compile()
#' dl=stan_format(R=data$recruits,S=data$spawners,t=data$broodyear)
stan_format<-function(s,r,t,K=2,smax_priors=NA){
  if(is.na(smax_priors)==T){
    dl=list(S=s,
            R_S=log(r/s),
            L=max(t)-min(t)+1,
            ii=t-min(t)+1,
            pSmax_mean=0.5*max(s),
            pSmax_sig=max(S),
            K=K,
            alpha_dirichlet=matrix(rep(1,K*2),ncol=K,nrow=K)
    )  
  }
  if(is.na(smax_priors)==F){
    dl=list(S=s,
            R_S=log(r/s),
            L=max(t)-min(t)+1,
            ii=t-min(t)+1,
            pSmax_mean=smax_priors[1],
            pSmax_sig=smax_priors[2],
            K=K,
            alpha_dirichlet=matrix(rep(1,K*2),ncol=K,nrow=K)
    )  
  }
  
} 


#' Simple Ricker model estimated with stan implemented in cmdstanr
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param lambertW Logical, indicating if lambertW functions should be used in stan code, requires git installation from stan
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
ricker_stan <- function(data,  AC=FALSE, control = stancontrol(), mod=NULL, warmup=300, chains = 6, iter = 1000,lambertW=FALSE,...) {
  
  
  if(is.null(mod)){
    sm <- compile_code(type='static',ac=AC,par='n',lambertW=lambertW)
  }else{
    sm <-mod
  }
  
  
  if(AC){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=seq_len(nrow(data)),
                R_S =data$logRS,
                S=data$S)
  }else{
    datm = list(N=nrow(data),
                R_S =data$logRS,
                S=data$S)
    
    
  }
  
  fit<-rstan::sampling(sm, data=datm,
                       control = control, warmup = warmup, 
                       chains = chains, iter = iter,verbose=FALSE)
  #fit <- rstan::stan(model_code = sm, 
  #                      data = datm,
  #                      control = control, warmup = warmup, chains = chains, iter = iter)
  #
  
  mc <- rstan::extract(fit, 
                       inc_warmup=FALSE, permuted=FALSE)
  
  
  aa <- rstan::summary(fit)
  
  ans<-list(alpha=aa$summary["log_a","50%"],
            beta=aa$summary["b","50%"],
            Smax=aa$summary["S_max","50%"],
            sigobs=aa$summary["sigma","50%"], 
            #Smsy=aa$summary["S_msy","50%"],
            #umsy=aa$summary["U_msy","50%"],
            stanfit=fit, 
            mcmcsummary=aa$summary,
            c_mcmcsummary=aa$c_summary, 
            samples=mc ) 
  
  if(lambertW){
    ans$Smsy=aa$summary["S_msy","50%"]
    ans$umsy=aa$summary["U_msy","50%"]
  }
  
  return(ans )
  
}