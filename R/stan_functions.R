# Functions to run the stan stock recruitment  models
#===============================================





#' Bayesian Ricker model with Stan
#'
#' @export
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
rickerstan <- function(data,...){
  standata <- list(TT=length(data$S),
  				R_S = data$logRS,
  				S = data$S)
  	
  out <- rstan::sampling(stanmodels$ricker_linear, data = standata, ...)
  return(out)
}



#==============================================================


# END
#***********************************************************************************
