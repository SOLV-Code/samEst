# Functions to run the stan stock recruitment  models
#===============================================





#' Bayesian linear regression with Stan
#'
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @export
#' 
rickerstan <- function(data,...) {
  standata <- list(TT=length(data$S),
  				R_S = data$logRS,
  				S = data$S)
  	
  out <- rstan::sampling(stanmodels$ricker_linear, data = standata, ...)
  return(out)
}



#==============================================================
#Dummy function

#' stanRoxygen commands
#'
#' This is a dummy function to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon package install.
#'
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @useDynLib samEst, .registration = TRUE
#'
#' 
#' 
#'
dummystan <- function(){
  return(NULL)
}

# END
#***********************************************************************************
