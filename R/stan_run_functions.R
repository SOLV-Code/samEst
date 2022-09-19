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
#' @importFrom rstan sampling 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerstan(data=harck)
#' 
rickerstan <- function(data,  AC=FALSE, control = stancontrol(),  chains = 6, iter = 1000,...) {

  sm <- sr_mod(type='static',ac=AC,par='n',loglik=FALSE,modelcode=TRUE)
  
  fit <- rstan::stan(model_code = sm, 
                        data = list(N=nrow(data),
                                    R_S =data$logRS,
                                    S=data$S),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  print(fit)
  aa<-rstan::summary(fit)
  


  return(list(alpha=aa$summary["log_a","50%"],
   beta=aa$summary["b","50%"],
   sigobs=aa$summary["sigma_e","50%"], 
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=samples ) )

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



