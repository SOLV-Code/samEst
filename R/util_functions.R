# util functions
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Find order of magnitude
#'
#' @param x a number. 
#' 
#' @returns order of magnitude
#' 
#' 
log10_ceiling <- function(x) {
    10^(ceiling(log10(x)))
}

log10_ceiling(max(data$S))



# END
#***********************************************************************************
