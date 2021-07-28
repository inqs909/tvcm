#' summary.tvcm: Plots time-varying coefficient models
#' 
#' This function allows you to plot the results from any of time-varying coefficient functions.
#'
#' @param object An object of class "tvcm"
#' @param bounds If TRUE, the bootstrap percentiles will be plotted
#' @param deriv If TRUE, the derivatives will be plotted
#' @param ... further arguments; not in use
#'
#'
#' @return Plots
#'  
#' 
#' @export
#'
#' @examples
#' 
#' fit <- tvcm(formula = Y~x1+x2, data = normal_data, time = time,
#'             id = id, se = TRUE, nboot = 100)
#'             
#'            
#' summary.tvcm(fit, bounds = TRUE, deriv = TRUE)
#' 
#' 
summary.tvcm <- function(object, bounds = FALSE, deriv = FALSE, ...){
  tp <- object$time_points
  est <- object$estimates$estimates
  nc <- rownames(est)
  ll <- length(nc)
  if (isFALSE(bounds)){
    for (i in 1:ll){
      plot(x = tp, y = est[i,],
           main = paste("Time-Varying for ", nc[i], sep = ""),
           xlab = "Time",
           ylab = "Estimate",
           type = "l"
      )
    }
  } else {
    if (is.null(object$bootstrap_results)){stop("Bootstrap Errors not computed")}
    lower <- object$bootstrap_results$boot_lower$est
    upper <- object$bootstrap_results$boot_upper$est
    for (i in 1:ll){
      minp <- min(lower[i,]) - 2 * abs(min(dist(lower[i,])))
      maxp <- max(upper[i,]) + 2 * abs(min(dist(upper[i,])))
      plot(x = tp, y = est[i,],
           main = paste("Time-Varying for ", nc[i], sep = ""),
           xlab = "Time",
           ylab = "Estimate",
           type = "l",
           ylim = c(minp, maxp)
      )
      graphics::lines(x = tp, y = lower[i,],type = "s")
      graphics::lines(x = tp, y = upper[i,],type = "s")
    }
  }
  if (isTRUE(deriv)){
    der <- object$estimates$deriv
    if (isFALSE(bounds)){
      for (i in 1:ll){
        plot(x = tp, y = der[i,],
             main = paste("Time-Varying Derivative for ", nc[i], sep = ""),
             xlab = "Time",
             ylab = "Estimate",
             type = "l"
        )
      }
    } else {
      if (is.null(object$bootstrap_results)){stop("Bootstrap Errors not computed")}
      lowerd <- object$bootstrap_results$boot_lower$deriv
      upperd <- object$bootstrap_results$boot_upper$deriv

      for (i in 1:ll){
        minpd <- min(lowerd[i,]) - 2 * abs(min(dist(lowerd[i,])))
        maxpd <- max(upperd[i,]) + 2 * abs(min(dist(upperd[i,])))
        plot(x = tp, y = der[i,],
             main = paste("Time-Varying Derivative for ", nc[i], sep = ""),
             xlab = "Time",
             ylab = "Estimate",
             type = "l",
             ylim = c(minpd, maxpd)
        )
        graphics::lines(x = tp, y = lowerd[i,],type = "s")
        graphics::lines(x = tp, y = upperd[i,],type = "s")
      }
    }
  }
}

