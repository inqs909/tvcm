#' kernel_cpp
#' 
#' Obtains id of kernel function for cpp
#'
#' @param type name of kernel
#'
#' @export
#'

kernel_cpp <- function(type){
  if (type %in% c("normal","epanechnikov","triweight","quartic")){
    if (type == "epanechnikov"){
      post = 1
    } else if (type == "normal") {
      post = 2
    } else if (type == "triweight") {
      post = 3
    } else if (type == "quartic") {
      post = 4
    }
  } else {
    stop("type must be normal, triweight, quartic, or epanechnikov")
  }
  return(as.integer(post))
}