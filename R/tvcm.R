tvcm <- function(formula, data, time, id = NULL,
                       method = "local-linear",
                       local_params = NULL,
                       p_splines_params = NULL,
                       bayesian_p_splines_params = NULL,
                       se = FALSE, alpha = 0.05, nboot = 1000){
  if (method == "local-linear"){
    if (is.null(local_params)){
      results <- tvcm_local(formula = formula, data = data, time = time, 
                 id = id, se = se, alpha = alpha, nboot = nboot)
      return(results)
    } else {
      params$ngrids <- ifelse(is.null(local_params$ngrids), 200, local_params$ngrids)
      params$gp <- ifelse(is.null(local_params$grid_points), NULL, local_params$grid_points)
      params$bw <- ifelse(is.null(local_params$bandwidth), NULL, local_params$bandwidth)
      params$kernel <- ifelse(is.null(local_params$kernel), "epanechnikov", local_params$kernel)
      results <- tvcm_local(formula = formula, data = data, time = time, id = id,
                            ngrids = params$ngrids, grid_points = params$gp,
                            bandwidth = params$bw, kernel = params$kernel,
                            se = se, alpha = alpha, nboot = nboot)
      return(results)
    }

  } 
}