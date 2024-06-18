


#' Calculate effective degrees of freedom code from 
#' Thorson JT. Measuring complexity for hierarchical models using effective 
#' degrees of freedom. Ecology. 2024 Jun 11:e4327. doi: 10.1002/ecy.4327. 
#' Epub ahead of print. PMID: 38859712.
#'
#' @param obj tmb_obj output from samEst TMB function
#' @param opt model output from samEst TMB function
#' @param nonvariance_fixed_effects names of the parameters that are the non varinace 
#' fixed effects in models with random effects
#' @param prediction_name predicted variable name
#' @param data_name observed y value
#' @param delta see thorson paper
#' @param refit refit for all data points
#' @param what what to output
#' 
#' @export
#' 
#' @returns Sgen likelihood 
#' 
#' 
#' 
calculate_EDF <- function( obj,
          opt,
          nonvariance_fixed_effects = "",
          prediction_name = "yhat_t",
          data_name = "obs_logRS",
          delta = 0.01,
          show_progress = FALSE,
          refit = c("random", "full")[1],
          what = c("EDF","grad_i")[1] ){

  # Extract stuff
  pred = obj$report()[[prediction_name]]
  obs = obj$env$data[[data_name]]

  # Calculate
  # get grid of runs
  if( is.array(pred) ){
    Dim = dim(pred)
  }else{
    Dim = length(pred)
  }
  grid_iz = as.matrix(expand.grid(lapply(Dim,FUN=seq_len)))
  pred_i = pred[grid_iz]

  # Optional progress bar
  if(isTRUE(show_progress)){
    require(progress)
    pb = progress_bar$new(total = nrow(grid_iz))
  }

  # Loop through data
  prednew_i = grad_i = rep(NA,nrow(grid_iz))
  for( i in seq_along(grad_i) ){
    if(isTRUE(show_progress)) pb$tick()
    # Update values
    obj$env$data[[data_name]][grid_iz[i,,drop=FALSE]] = obj$env$data[[data_name]][grid_iz[i,,drop=FALSE]] + delta
    # Re-run
    if( refit=="full"){
      opt_t = nlminb( start=opt$par, obj=obj$fn, grad=obj$gr )
    }else if( refit=="random" ){
      obj$fn(opt$par)
    }else stop("Check `refit`")
    # Calculate gradient
    prednew = obj$report()[[prediction_name]]
    prednew_i[i] = prednew[grid_iz[i,,drop=FALSE]]
    grad_i[i] = (prednew_i[i] - pred_i[i]) / delta
    # Reset values
    obj$env$data[[data_name]] = obs
  }

  # Compute variance params
  extra_params = names(opt$par)
  if( refit=="full" ){
    matches = which( extra_params %in% nonvariance_fixed_effects )
    if(length(matches)>0) extra_params = extra_params[-matches]
  }

  # Compute degrees of freedom
  EDF = sum(grad_i) + length(extra_params) # Plus add any variance parameters + 1
  if(what=="EDF") return( EDF )
  if(what=="grad_i") return( grad_i )
  if(what=="all") return( list("EDF"=EDF, "grad_i"=grad_i, "prednew_i"=prednew_i, "pred_i"=pred_i) )
}
