#Function to calculate reference points once SR parameters are computed. 






#sgen functions from samSim



#' compute Sgen likelihood function
#'
#' @param S Spawner numbers set to the interval between 0 and Smsy in sGenSolver.
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' 
#' 
#' @returns Sgen likelihood 
#' 
#' 
#' 
Sgencompute <- function(S, a,b, Smsy ) {
	#modified from samsim sGenOptimSmsyum
  
  prt <- S * exp(a - b * S)
  epsilon <- log(Smsy) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, 0.01, log = T))
  return(list(nSS = sum(nLogLike)))
  
}





#' compute Sgen estimate
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' 
#' 
#' @returns Sgen estimate 
#' 
#' 
#' 
sGenSolver <- function(a,b, Smsy) {
  #gives the min Ricker log-likelihood
  if(a>0){
    fnSGen <- function(S, a, b, Smsy) -1.0 * Sgencompute(S, a, b, Smsy)$nSS
    fit <- optimize(f = fnSGen, interval = c(0, ((a / b) * (0.5 - 0.07 * a))),
                 a=a,b=b, Smsy = Smsy)
  }else{
    fit <-list(minimum=NA)
  }

  return(list(fit = fit$minimum))
}






#' compute Smsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' 
#' 
smsySolver <- function(a,b) {
  #gives the min Ricker log-likelihood
  Smsy <- (1 - gsl::lambert_W0(exp(1 - a))) /b

  return(Smsy)
}



#' compute Umsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' 
#' 
umsySolver <- function(a,b) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - a)))

  return(umsy)
}
