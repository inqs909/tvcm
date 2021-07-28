#' tvcm: Time-Varying Coefficient Models
#' 
#' This function estimates the time-varying coefficient model at different time points.
#'
#'
#' @param formula A formula class object. Provides information about the names of the response and predictor variables.
#' @param data A data frame containing the variables needed for the function.
#' @param time The name of the variable indicating the time points for each observation. 
#' @param id The name of the variable to group observations. Necessary when dealing with repeated measurements such as longitudinal data.
#' @param ngrids When specified, a vector of size ngrid is created for the varying-coefficient values. The default value is 200. The vector is creates using the maximum and minimum from the provided time points.
#' @param grid_points A vector indicating the grid points to estimate the varying coefficient values. When specified, ngrid is ignored.
#' @param bandwidth A numeric value indicating the bandwidth. Default is the tenth of the range of the time points. 
#' @param kernel The name of the kernel function used to estimate the varying-coefficient values. Default is Epanechnikov.
#' @param se If set TRUE, a bootstrap method is applied to estimate the standard errors and percentiles. 
#' @param alpha A value indicating the significance level for the percentiles.
#' @param nboot A number indicating how many boot samples to construct.
#'
#' @return tvcm returns a list containing the estimated varying coefficients
#' \itemize{
#' \item estimates: A list containing two matrices:
#' \itemize{
#' \item est: a matrix for the estimates, each row represents the estimate of the varying coefficient
#' \item deriv: a matrix for the first derivative, each row represents the derivative of the varying coefficient
#' }
#' \item bootstrap_results: If specified, a list containing the elements below
#' \itemize{
#' \item boot_se: A list containing two matrices for the standard errors:
#' \itemize{
#' \item est: a matrix for the estimates, each row represents the estimate of the varying coefficient
#' \item deriv: a matrix for the first derivative, each row represents the derivative of the varying coefficient
#' }
#' \item boot_lower: A list containing two matrices for the lower percentiles:
#' \itemize{
#' \item est: a matrix for the estimates, each row represents the estimate of the varying coefficient
#' \item deriv: a matrix for the first derivative, each row represents the derivative of the varying coefficient
#' }
#' \item boot_upper:A list containing two matrices for the upper percentiles:
#'  \itemize{
#' \item est: a matrix for the estimates, each row represents the estimate of the varying coefficient
#' \item deriv: a matrix for the first derivative, each row represents the derivative of the varying coefficient
#' }
#' \item boot_samples: An array containing the bootstrap samples
#' }
#' \item time_points: The time points used for the estimates
#' }
#' 
#' @export 
#' 
#' @author Isaac Quintanilla Salinas
#'
#' @examples 
#' 
#' tvcm(formula = Y~x1+x2, data = normal_data, time = time, 
#'      id = id, se = TRUE, nboot = 100)
#' 
tvcm <- function(formula, data, time, id = NULL,
                 ngrids = 200, grid_points = NULL,
                 bandwidth = NULL, kernel = "epanechnikov",
                 se = FALSE, alpha = 0.05, nboot = 1000){
  # Kernel Information ####
  kernel_id <- kernel_cpp(kernel) # ID for cpp
  
  # Model Information ####
  mf <- model.frame(formula, data)
  mr <- as.vector(model.response(mf))
  mm <- as.matrix(model.matrix(formula, mf))
  mt <- as.vector(data$time)
  dd <- dim(mm)[2]*2
  
  lvcm <- dim(mm)[2]
  ipp <- cbind(diag(rep(1,lvcm)), diag(rep(0,lvcm)))
  ipp_deriv <- cbind(diag(rep(0,lvcm)),diag(rep(1,lvcm)))
  
  # Setting Up Gridpoints ####
  if (is.null(grid_points)){
    time_range <- range(mt)
    gp <- seq(time_range[1], time_range[2], length.out = ngrids)
    gp_labels <- as.character(round(gp,2))
    ll <- ngrids
  } else {
    time_range <- range(mt)
    gp <- grid_points
    gp_labels <- as.character(round(gp,2))
    ll <- length(gp)
  }
  
  # Setting Up Bandwidth ####
  if (is.null(bandwidth)){
    hh <- (time_range[2] - time_range[1])/10
  } else {
    hh <- bandwidth
  }
  # Estimating VCMs####
  if (is.null(data$id)) {
    
    rvcm<-sapply(gp, vcm_wls,
           x = mm, y = mr, time = mt, h= hh, type = kernel_id,
           simplify = "matrix")
    
  } else {
    mid <- as.vector(data$id)
    rvcm<-sapply(gp, tvcm_wls,
                 x = mm, y = mr, time = mt, id = mid, h= hh, type = kernel_id,
                 simplify = "matrix")
    
  }
  # #Estimating Sigmas####
  # if (isFALSE(sigma)){
  #   sigmas <- NULL
  # } else {
  #   tll <- length(mt)
  #   ktt0 <- matrix(nrow = tll , ncol = ll)
  #   for (i in 1:ll){
  #     ktt0[, i] <- kernel_vec(mt, gp[i], hh, kernel_id)
  #   }
  #   if (is.null(data$id)){
  #     smt<-unique(mt)
  #     sll <- length(smt)
  #     srvcm <- matrix(nrow = sll, ncol = dd)
  #     for (i in 1:sll){
  #       rvcm[, i] <- vcm_wls(mm, mr, mt, smt[i], hh, kernel_id)
  #     }
  #     stid <- vector(length = tll)
  #     for (i in 1:tll){
  #       stid[i] <- which(mt[i]==smt)
  #     }
  #     sigmas <- vector(length = ll)
  #     res2 <- vector(length = tll)
  #     for (i in 1:tll){
  #       res2[i] <- (mr[i] - mm[i,] %*% ipp %*% rvcm[stid[i],]) ^ 2
  #     }
  #     for (i in 1:ll){
  #       sigmas[i]<-res2 %*% ktt0[, i] / sum(ktt0[, i])
  #     }
  #   } else {
  #     smt<-unique(mt)
  #     sll <- length(smt)
  #     mid <- as.vector(data$id)
  #     srvcm <- matrix(nrow = sll, ncol = dd)
  #     for (i in 1:sll){
  #       srvcm[i, ] <- tvcm_wls(mm, mr, mt, mid, smt[i], hh, kernel_id)
  #     }
  #     stid <- vector(length = tll)
  #     for (i in 1:tll){
  #       stid[i] <- which(mt[i]==smt)
  #     }
  #     sigmas <- vector(length = ll)
  #     res2 <- vector(length = tll)
  #     for (i in 1:tll){
  #       res2[i] <- (mr[i] - mm[i,] %*% ipp %*% rvcm[stid[i],]) ^ 2
  #     }
  #     for (i in 1:ll){
  #       sigmas[i]<-res2 %*% ktt0[, i] / sum(ktt0[, i])
  #     }
  #   }
  # }
  # 
  # Obtaining Bootstrap Standard Errors ####
  
  if (se) {
    bdim <- dim(rvcm)
    if (is.null(data$id)) {
      boot_max <- length(mr)
      boots_array <- array(dim = c(bdim, nboot))
      for (i in 1:nboot){
        boot_ids <- sample(1:boot_max, boot_max, replace = TRUE)
        boot_data <- c()
        for (ii in boot_ids){
          boot_sub <- data[ii,]
          boot_data <- rbind(boot_data, boot_sub)
        }
        boot_df <- as.data.frame(boot_data)
        bmf <- model.frame(formula, boot_df)
        bmr <- as.vector(model.response(bmf))
        bmm <- as.matrix(model.matrix(formula, bmf))
        bmt <- as.vector(boot_df$time)
        
        brvcm <- matrix(nrow = ll, ncol = dd)

        brvcm<-sapply(gp, vcm_wls,
                      x = bmm, y = bmr, time = bmt, h= hh, type = kernel_id,
                      simplify = "matrix")
        
        boots_array[,,i] <- brvcm
      }
      boot_se <- apply(boots_array, c(1,2), sd)
      boot_ll <- apply(boots_array, c(1,2), quantile, probs = alpha/2)
      boot_ul <- apply(boots_array, c(1,2), quantile, probs = 1-alpha/2)
      
      boot_res <- list(boot_se = boot_se, boot_lower = boot_ll,
                       boot_upper = boot_ul, significance_level = alpha,
                       boot_samples = boots_array)
    } else {
      boot_max <- max(data$id)
      boots_array <- array(dim = c(bdim, nboot))
      for (i in 1:nboot){
        boot_ids <- sample(1:boot_max, boot_max, replace = TRUE)
        boot_data <- c()
        for (ii in boot_ids){
          boot_sub <- subset(data, data$id == ii)
          boot_data <- rbind(boot_data, boot_sub)
        }
        boot_df <- as.data.frame(boot_data)
        bmf <- model.frame(formula, boot_df)
        bmr <- as.vector(model.response(bmf))
        bmm <- as.matrix(model.matrix(formula, bmf))
        bmt <- as.vector(boot_df$time)
        bmid <- as.vector(boot_df$id)
        
        brvcm<-sapply(gp, tvcm_wls,
                      x = bmm, y = bmr, time = bmt, id = bmid, h= hh, type = kernel_id,
                      simplify = "matrix")
        
        boots_array[,,i] <- brvcm
      }
      boot_se <- apply(boots_array, c(1,2), sd)
      boot_ll <- apply(boots_array, c(1,2), quantile, probs = alpha/2)
      boot_ul <- apply(boots_array, c(1,2), quantile, probs = 1-alpha/2)
      
      
      boot_se_est <- ipp %*% boot_se
      boot_se_deriv <- ipp_deriv %*% boot_se
      
      boot_ll_est <- ipp %*% boot_ll
      boot_ul_est <- ipp %*% boot_ul
      
      boot_ll_deriv <- ipp_deriv %*% boot_ll
      boot_ul_deriv <- ipp_deriv %*% boot_ul
      
      rownames(boot_se_est) <-  colnames(mm)
      rownames(boot_se_deriv) <- colnames(mm)
      
      rownames(boot_ll_est) <-  colnames(mm)
      rownames(boot_ll_deriv) <- colnames(mm)
      
      rownames(boot_ul_est) <-  colnames(mm)
      rownames(boot_ul_deriv) <- colnames(mm)
      
      colnames(boot_se_est) <-  gp_labels
      colnames(boot_se_deriv) <- gp_labels
      
      colnames(boot_ll_est) <-  gp_labels
      colnames(boot_ll_deriv) <- gp_labels
      
      colnames(boot_ul_est) <-  gp_labels
      colnames(boot_ul_deriv) <- gp_labels
      
      
      boot_res <- list(boot_se = list(est = boot_se_est, deriv = boot_se_deriv), boot_lower = list(est = boot_ll_est, deriv = boot_ll_deriv),
                       boot_upper = list(est = boot_ul_est, deriv = boot_ul_deriv), significance_level = alpha,
                       boot_samples = boots_array)
      
    }
  } else {
    boot_res <- NULL
  }
  
  rvcm_est <- ipp %*% rvcm
  rvcm_deriv <- ipp_deriv %*% rvcm 
  
  rownames(rvcm_est) <-  colnames(mm)
  rownames(rvcm_deriv) <- colnames(mm)
  
  colnames(rvcm_est) <- gp_labels
  colnames(rvcm_deriv) <- gp_labels
  
  
  results <- list(estimates = list(estimates = rvcm_est, deriv = rvcm_deriv), bootstrap_results = boot_res, time_points = gp)
  class(results) <- "tvcm"
  
  return(results)
  
}