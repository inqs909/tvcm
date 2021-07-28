#' summary.tvcm.cv
#' 
#' Prints the recommended bandwidth.
#'
#' @param object An object of class "tvcm.cv"
#' @param ... further arguments; not in use
#'
#' @export
#'
#' 
summary.tvcm.cv <- function(object, ...){
  print(paste("Recommended Bandwidth: ",object$recommended_bandwidth, sep = ""))
  print(paste("Negative Log-Likelihood Value: ",object$value, sep = ""))
}