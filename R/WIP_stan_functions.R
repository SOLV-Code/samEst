#' stan compile - compiles stan versions of static, autocorrelated, random walk, and regime shift models
#'
#' @param type Options for:
#' * all (default) - compiles 8 model sets
#' * static - compiles the classic (static) Ricker model
#' * ac - compiles the static Ricker model with 1-y residual autocorrelation
#' * tv - 'time-varying' model, must also specify the parameter ('par') - see below
#' * regime - regime shift (hidden markov) model, must also specify the parameter ('par') - see below
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
  if(file.exists(file.path(cmdstanr::cmdstan_path(),'samEst'))==FALSE){
    dir.create(file.path(cmdstanr::cmdstan_path(),'samEst'))
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[1]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[1])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[2]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[2])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[3]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[3])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[4]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[4])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[5]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[5])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[6]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[6])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[7]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[7])
    cmdstanr::write_stan_file(code=readLines(file.path('./src/stan/',list.files('./src/stan/','.stan'))[8]),dir=file.path(cmdstanr::cmdstan_path(),'samEst'),basename=list.files('./src/stan/','.stan')[8])
    }
  
  if(type=='all'){
   m1=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[1]),pedantic=F)
   m2=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[2]),pedantic=F)
   m3=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[3]),pedantic=F)
   m4=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[4]),pedantic=F)
   m5=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[5]),pedantic=F)
   m6=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[6]),pedantic=F)
   m7=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[7]),pedantic=F)
   m8=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[8]),pedantic=F)
   
  }else if(type=='static'){
    m1=cmdstanr::cmdstan_model(exe_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.exe')[1]),pedantic=F)
  }else if(type=='ac'){
    m2=cmdstanr::cmdstan_model(exe_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.exe')[2]),pedantic=F)
  }else if(type=='tv'&par=='prod'){
    m3=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[3]),pedantic=F)
  }else if(type=='tv'&par=='cap'){
    m4=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[4]),pedantic=F)
  }else if(type=='tv'&par=='both'){
    m5=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[5]),pedantic=F)
  }else if(type=='regime'&par=='prod'){
    m6=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[6]),pedantic=F)
  }else if(type=='regime'&par=='cap'){
    m7=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[7]),pedantic=F)
  }else if(type=='regime'&par=='both'){
    m8=cmdstanr::cmdstan_model(stan_file=file.path(cmdstanr::cmdstan_path(),'samEst',list.files(file.path(cmdstanr::cmdstan_path(),'samEst'),'.stan')[8]),pedantic=F)
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
            pSmax_sig=max(s),
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
ricker_stan<- function(s,r,t,K=2,smax_priors=NA,  control = stancontrol(), type=c('static','ac','tv','regime'),par=c('prod','cap','both'), warmup=300, chains = 6, iter = 1000){
  
  dl=stan_format(s=s,r=r,t=t,K=K,smax_priors=smax_priors)
  
  if(type=='static'){
    if(exists('m1')==FALSE){
      stan_compile(type='static')
    }else{
      
      fit=m1$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }else if(type=='ac'){
    if(exists('m2')==FALSE){
      stan_compile(type='ac')
    }else{
      fit=m2$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='tv'&par=='prod'){
    if(exists('m3')==FALSE){
      stan_compile(type='tv',par='prod')
    }else{
      fit=m3$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='tv'&par=='cap'){
    if(exists('m4')==FALSE){
      stan_compile(type='tv',par='cap')
    }else{
      fit=m4$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='tv'&par=='both'){
    if(exists('m5')==FALSE){
      stan_compile(type='tv',par='both')
    }else{
      fit=m5$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='regime'&par=='prod'){
    if(exists('m6')==FALSE){
      stan_compile(type='regime',par='prod')
    }else{
      fit=m6$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='regime'&par=='cap'){
    if(exists('m7')==FALSE){
      stan_compile(type='regime',par='cap')
    }else{
      fit=m7$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  if(type=='regime'&par=='both'){
    if(exists('m8')==FALSE){
      stan_compile(type='regime',par='both')
    }else{
      fit=m8$sample(data=dl,
                    chains=chains,
                    iter_warmup=warump,
                    iter_sampling=iter,
                    refresh=0,
                    adapt_delta=0.999,
                    max_treedepth=20)
    }
  }
  
  return(fit)
}

ref_points<-function(fit,SGen=FALSE,plot=c('all','alpha','Smax','Smsy','Umsy','Sgen')){
  
  if(SGen==TRUE){
    df=data.frame(par=c('log_alpha','Smax','Smsy','Umsy','Sgen'),median=NA,mean=NA,sd=NA,l95=NA,u95=NA,l80=NA,u80=NA)
    
    draws=fit$draws(variables=c('log_a','b','Smax','Smsy','Umsy'),format='draws_matrix')
    df[1,2]=median(draws[,1]);df[1,3]=mean(draws[,1]);df[1,4]=sd(draws[,1]);df[1,5]=quantile(draws[,1],0.025);df[1,6]=quantile(draws[,1],0.975);df[1,7]=quantile(draws[,1],0.1);df[1,8]=quantile(draws[,1],0.9)
    df[2,2]=median(draws[,3]);df[2,3]=mean(draws[,3]);df[2,4]=sd(draws[,3]);df[2,5]=quantile(draws[,3],0.025);df[2,6]=quantile(draws[,3],0.975);df[2,7]=quantile(draws[,3],0.1);df[2,8]=quantile(draws[,3],0.9)
    df[3,2]=median(draws[,4]);df[3,4]=mean(draws[,4]);df[3,4]=sd(draws[,4]);df[3,5]=quantile(draws[,4],0.025);df[3,6]=quantile(draws[,4],0.975);df[3,7]=quantile(draws[,4],0.1);df[3,8]=quantile(draws[,4],0.9)
    df[4,2]=median(draws[,4]);df[4,3]=mean(draws[,4]);df[4,4]=sd(draws[,4]);df[4,5]=quantile(draws[,4],0.025);df[4,6]=quantile(draws[,4],0.975);df[4,7]=quantile(draws[,4],0.1);df[4,8]=quantile(draws[,4],0.9)
    df[5,2]=median(draws[,5]);df[5,3]=mean(draws[,5]);df[5,4]=sd(draws[,5]);df[5,5]=quantile(draws[,5],0.025);df[5,6]=quantile(draws[,5],0.975);df[5,7]=quantile(draws[,5],0.1);df[5,8]=quantile(draws[,5],0.9)
  
    sgen=apply(draws,1,samEst::sGenCalc,a=draws[,1],b=draws[,2],Smsy=draws[,4])
    
    df[1,2]=median(sgen);df[1,3]=mean(sgen);df[1,4]=sd(sgen);df[1,5]=quantile(sgen,0.025);df[1,6]=quantile(sgen,0.975);df[1,7]=quantile(sgen,0.1);df[1,8]=quantile(sgen,0.9)
    
  }
  if(SGen==FALSE){
    df=data.frame(par=c('log_alpha','Smax','Smsy','Umsy'),median=NA,mean=NA,sd=NA,l95=NA,u95=NA,l80=NA,u80=NA)
  }
      
  
  
  
return(df)
}
