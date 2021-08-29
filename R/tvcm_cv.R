tvcm_cv <- function(formula, data, time, id = NULL,
                    method = "local-linear",
                    local_params = NULL,
                    p_splines_params = NULL,
                    bayesian_p_splines_params = NULL){
  if (method == "local-linear"){
    if (is.null(local_params)){
      results <- tvcm_cv_local(formula = formula, data = data,
                               time = time, id = id)
      return(results)
    } else {
      params$bw <- ifelse(is.null(local_params$bandwidth), NULL, local_params$bandwidth)
      params$kernel <- ifelse(is.null(local_params$kernel), "epanechnikov", local_params$kernel)
      results <- tvcm_cv_local(formula = formula, data = data, time = time, id = id,
                               bandwidth = params$bw, kernel = params$kernel)
      return(results)
    }
  }
}