#' A simulated dataset of 200 observations with each containing 25 correlated normal response variables.
#'
#'
#' @format A data frame with 5000 rows and 5 variables:
#' \describe{
#'   \item{id}{A variable grouping the repeated measurements}
#'   \item{time}{The time for each repeated measurement}
#'   \item{x1}{First predictor variable}
#'   \item{x2}{Second predictor variable}
#'   \item{Y}{Correlated normal response variable}
#'  }
#'  
#'
"normal_data"

#' A simulated dataset of 200 observations with each containing 25 correlated binary response variables.
#'
#'
#' @format A data frame with 5000 rows and 5 variables:
#' \describe{
#'   \item{id}{A variable grouping the repeated measurements}
#'   \item{time}{The time for each repeated measurement}
#'   \item{x1}{First predictor variable}
#'   \item{x2}{Second predictor variable}
#'   \item{Y}{Correlated binary response variable}
#'  }
#'  
#'  
"binary_data"


#' A simulated dataset of 200 observations with each containing 25 correlated Poisson response variables.
#'
#'
#' @format A data frame with 5000 rows and 5 variables:
#' \describe{
#'   \item{id}{A variable grouping the repeated measurements}
#'   \item{time}{The time for each repeated measurement}
#'   \item{x1}{First predictor variable}
#'   \item{x2}{Second predictor variable}
#'   \item{Y}{Correlated Poisson response variable}
#'  }
#'  
#'  
"pois_data"