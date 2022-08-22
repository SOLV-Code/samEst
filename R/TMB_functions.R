# Functions to run the TMB stock recruitment  models
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Simple Ricker model fir with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' 
#' 
#' 
#' @export
rickerTMB <- function(data,  silent = FALSE, control = TMBcontrol()) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    obs_S = data$S,
    obs_logRS = data$logRS
  )

  initlm<-lm(logRS~S, data=data)

  tmb_params <- list(
    alpha   = initlm$coefficients[[1]],
    logbeta = ifelse(initlm$coefficients[[2]]<0,1e-08,log(-initlm$coefficients[[2]])),
    logsigobs = log(1)
  )

  #to be implemented
  tmb_map <- list()
  tmb_random <- NULL

  #===================================
  # TMB fit
  #===================================

  tmb_obj <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params, map = tmb_map,
      random = tmb_random, DLL = "Ricker_simple", silent = silent)
  
  tmb_opt <- stats::nlminb(
    start = tmb_obj$par, objective = tmb_obj$fn, gradient = tmb_obj$gr,
    control = control)
  
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)

  #todo add alpha, beta and sigma parameter esitimates

  structure(list(
    alpha    = tmb_obj$report()$alpha,
    beta     = tmb_obj$report()$beta,
    sig      = tmb_obj$report()$sigobs,
    model      = tmb_opt,
    data       = data,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_random = tmb_random,
    tmb_obj    = tmb_obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = "Ricker_simple")

}






#' Optimization control options. Copied from sdmTMB
#'
#' Any arguments to pass to [stats::nlminb()].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @export
TMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}





#' get TMB convergence diagnostics taken from sdmTMB
#'
#' @param sd_report A TMB sd report object
#' @export
get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
        "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
          "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
          "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  invisible(list(final_grads = final_grads, bad_eig = bad_eig))
}



#==============================================================
#Dummy function

#' Roxygen commands
#'
#' This is a dummy function to hold the useDynLib roxygen tag.
#' This tag will populate the namespace with compiled c++ functions upon package install.
#'
#' @useDynLib Ricker_simple
#' @useDynLib Rickerkf
#'
dummy <- function(){
  return(NULL)
}

# END
#***********************************************************************************
