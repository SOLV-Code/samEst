#----------------------------------------------------
#Catarina put this functions together to be used in simulation
# 
#----------------------------------------------------




#' Simple Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
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
#' @importFrom reshape2 melt 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerstan(data=harck)
#' 
ricker_stan <- function(data,  AC=FALSE, control = stancontrol(), warmup=300,  chains = 6, iter = 1000,...) {
 
    sm <- sr_mod(type='static',ac=AC,par='n',loglik=FALSE,modelcode=TRUE)

  
  fit <- rstan::stan(model_code = sm, 
                        data = list(N=nrow(data),
                                    R_S =data$logRS,
                                    S=data$S),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  

    mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    mcmc<-reshape2::melt(mc, as.is=TRUE)
    
    aa<-rstan::summary(fit)
  


  return(list(alpha=aa$summary["log_a","50%"],
   beta=aa$summary["b","50%"],
   sigobs=aa$summary["sigma_e","50%"], 
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mcmc ) )

}





#' random walks Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
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
#' @importFrom reshape2 melt 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_rw_stan(data=harck)
#' 
ricker_rw_stan <- function(data, par=c('a','b','both'),  control = stancontrol(), warmup=300,  chains = 6, iter = 1000,...) {
  #par='b'
  sm <- sr_mod(type='rw',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)

  
  fit <- rstan::stan(model_code = sm, 
                        data = list(N=nrow(data),
                                    L=nrow(data),
                                    ii=seq_along(data$logRS),
                                    R_S =data$logRS,
                                    S=data$S),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  

    mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    mcmc<-reshape2::melt(mc, as.is=TRUE)
    
    aa<-rstan::summary(fit)
  


  return(list(alpha=aa$summary["log_a","50%"],
   beta=aa$summary["b","50%"],
   sigobs=aa$summary["sigma_e","50%"],
   siga=ifelse(par=="a"|par=="both",aa$summary["sigma_a","50%"]),
   sigb=ifelse(par=="b"|par=="both",aa$summary["sigma_a","50%"]), 
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mcmc ) )

}



#' random walks Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
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
#' @importFrom reshape2 melt 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_hmm_stan(data=harck)
#' 
ricker_hmm_stan <- function(data, par=c('a','b','both'), k_regime=2, 
  control = stancontrol(), warmup=300,  chains = 6, iter = 1000,...) {
  #par='both'
  sm <- sr_mod(type='hmm',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)

  
  fit <- rstan::stan(model_code = sm, 
                        data = list(N=nrow(data),
                                    R_S =data$logRS,
                                    S=data$S,
                                    K=k_regime,
                                    alpha_dirichlet=rep(1,k_regime)
                                    ),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  

    mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    mcmc<-reshape2::melt(mc, as.is=TRUE)
    
    aa<-rstan::summary(fit)
  
row.names(aa$summary)[grep("^b\\[",row.names(aa$summary))]

aa$summary[grep("zstar",row.names(aa$summary)),"50%"]
# these are meaningless need to calc for each mcmc draw check with Dan
#need regime and weighted alpha, beta, sigma  
   return(list(alpha=aa$summary[grep("log_a",row.names(aa$summary)),"50%"],
   beta=aa$summary[grep("^b\\[",row.names(aa$summary)),"50%"],
   sigobs=aa$summary["sigma","50%"],
   pi=aa$summary[grep("pi1",row.names(aa$summary)),"50%"],
   A=aa$summary[grep("A",row.names(aa$summary)),"50%"],
   probregime =matrix(aa$summary[grep("gamma\\[",row.names(aa$summary)),"50%"],ncol=k_regime, byrow=T),
   regime = aa$summary[grep("zstar",row.names(aa$summary)),"50%"]
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mcmc ) )

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



