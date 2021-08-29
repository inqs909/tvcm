#' tvcm_cv_local: Cross-Validation for Time-Varying Coefficient Model 
#'
#' This function conducts a leave-one-out cross-validation approach to determine the best bandwidth for the data.
#'
#' @param formula A formula class object. Provides information about the names of the response and predictor variables. 
#' @param data A data frame containing the variables needed for the function.
#' @param time The name of the variable indicating the time points for each observation. 
#' @param id The name of the variable to group observations. Necessary when dealing with repeated measurements such as longitudinal data.
#' @param bandwidths A vector containing the bandwidth values for the cross-validation. Default is tests 50 bandwidths from ranging from double the smallest time point difference to range of time points. 
#' @param kernel The name of the kernel function used to estimate the varying-coefficient values. Default is Epanechnikov.
#'
#' @return tvcm_cv_local returns a list containing the results from the cross-validation
#' \itemize{
#' \item recommended_bandwidth: The bandwidth with the smallest negative log-likelihood.
#' \item value: value of negative log-likelihood at recommended bandwidth
#' \item rss: The tested bandwidths and their corresponding negative log-likelihood value.
#' }
#' @export
#' 
#' @author Isaac Quintanilla Salinas
#' 
#' @examples 
#' cv_fit <- tvcm_cv_local(formula = Y~x1+x2, data = pois_data, time = time, id = id)
#' 
#' summary(cv_fit)
#' 
tvcm_cv_local <- function(formula, data, time, id = NULL,
                          bandwidths = NULL, kernel = "epanechnikov"){
  # Kernel Information
  kernel_id <- kernel_cpp(kernel) # ID for cpp

  if (is.null(bandwidths)){
    time_unique <- unique(data$time)
    time_range <- max(time_unique) - min(time_unique)
    time_min <- min(dist(time_unique))*2
    hh <- seq(time_min, time_range, length.out = 20)
  } else {
    hh <- bandwidths
  }
  
  res <- matrix(nrow = length(hh), ncol = 2)
  
  for (bb in 1:length(hh)){
    hhh <- hh[[bb]]
    # Setting Up Gridpoints
    if (is.null(data$id)) {
        id_max <- length(data$time)
        cv_rss <- vector(length = id_max)
        for (cv in 1:id_max){
          cvii <- data[cv,]
          cviimf <- model.frame(formula, cvii)
          cviimr <- as.vector(model.response(cviimf))
          cviimm <- as.matrix(model.matrix(formula, cviimf))
          cviimt <- as.vector(cvii$time)
          
          
          cv_df <- data[-cv,]
          cvmf <- model.frame(formula, cv_df)
          cvmr <- as.vector(model.response(cvmf))
          cvmm <- as.matrix(model.matrix(formula, cvmf))
          cvmt <- as.vector(cv_df$time)
          
          lvcm <- dim(cvmm)[2]
          ipp <- cbind(diag(rep(1,lvcm)), diag(rep(0,lvcm)))
          
  
          nvcm <- vcm_wls_local(cvmm, cvmr, cvmt,
                          as.numeric(cviimt), hhh, kernel_id)
          cv_rss[cv] <- (cviimr - cviimm %*% ipp %*% nvcm)^2
        }
        pre_results <-sum(cv_rss)      
        
      } else {
        cv_ids <- unique(data$id)
        id_max <- length(cv_ids)
        cv_rss <- vector(length = id_max)
        for (cv in 1:id_max){
          cvii <- subset(data, data$id==cv)
          cviimf <- model.frame(formula, cvii)
          cviimr <- as.vector(model.response(cviimf))
          cviimm <- as.matrix(model.matrix(formula, cviimf))
          cviimt <- as.vector(cvii$time)
          
          lvcm <- dim(cviimm)[2]
          ipp <- cbind(diag(rep(1,lvcm)), diag(rep(0,lvcm)))
          
          cv_df <- subset(data, data$id!=cv)
          cvmf <- model.frame(formula, cv_df)
          cvmr <- as.vector(model.response(cvmf))
          cvmm <- as.matrix(model.matrix(formula, cvmf))
          cvmt <- as.vector(cv_df$time)
          cvmid <- as.vector(cv_df$id)
          
          cvrvcm <- matrix(ncol = length(cvii$time), nrow = 2*lvcm)
          for (i in 1:length(cviimt)){
            cvrvcm[, i] <- tvcm_wls_local(cvmm, cvmr, cvmt, cvmid, cviimt[i], hhh, kernel_id)
          }
          cviirs <- vector(length = length(cviimt))
          for (i in 1:length(cviimt)){
            cviirs[i] <- (cviimr[i] - cviimm[i,] %*% ipp %*% cvrvcm[,i])^2
          }
          cv_rss[cv] <- sum(cviirs)
        }
        pre_results <- sum(cv_rss) 
      }
    res[bb,] <- c(hhh, pre_results)
  }  
  colnames(res) <- c("bandwidth", "neg loglik")
  b_id<-which.min(res[,2])
  band <- res[b_id,1]
  results <- list(recommended_bandwidth = band,
                  value = res[b_id,2],
                  neg_logliks = res,
                  method = "local-linear")
  class(results) <- "tvcm.cv"
  return(results)
}

