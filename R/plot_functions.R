#' sr_plot function
#'
#' This function generates a stock-recruitment (S-R) model for rstan based on model inputs.
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'tv' model, or a regime shift model 'regime'
#' @param df stock-recruitment dataset. With columns R (recruits) and S (spawners) and by (brood year)
#' @param form Either fit with Stan ('stan') or TMB ('tmb')
#' @param mod a fitted Stan or TMB model
#' @param pdf whether to create a pdf or not
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model
#' @export
#' @examples
#' sr_plot(type='static',df=df,form='stan',par='n',loglik=T)
sr_plot<- function(type=c('static','rw','hmm'),df,form=c('stan','tmb'),mod, pdf=0){
  if(pdf==1){
    
  }
  if(type=='static'){
    
    
  }
  
}