#' @export
summary.FitResult <- function(object) {
  if(any(class(object$mod) == "aghq")){
    aghq_summary <- summary(object$mod)
    aghq_output <- capture.output(aghq_summary)
    cur_index <- grep("Here are some moments and quantiles for the transformed parameter:", aghq_output)
    cat(paste(aghq_output[1:(cur_index - 1)], collapse = "\n"))
    
    if(class(aghq_summary) == "aghqsummary"){
      cat("\nHere are some moments and quantiles for the log precision: \n")
      # cur_index <- grep("median    97.5%", aghq_output)
      # cat(paste(aghq_output[cur_index], collapse = "\n"))
    }
    
    if(!is.null(aghq_summary$summarytable)){
      summary_table <- as.matrix(aghq_summary$summarytable)
      theta_names <- c()
      for (instance in object$instances) {
        theta_names <- c(theta_names, paste("theta(", toString(instance@smoothing_var), ")", sep = ""))
      }
      if((length(row.names(summary_table)) - length(theta_names)) >= 1 ){
        num_theta_family <- (length(row.names(summary_table)) - length(theta_names))
        theta_names <- c(theta_names, rep(paste("theta(", "family", ")", sep = ""),  num_theta_family))
      }
      row.names(summary_table) <- theta_names
      print(summary_table)
      cat("\n")
    }
  }
  # samps <- aghq::sample_marginal(object$mod, M = 3000)
  samps <- object$samps
  if(length(unlist(object$fixed_samp_indexes)) >= 1){
    fixed_samps <- samps$samps[unlist(object$fixed_samp_indexes), , drop = F]
    fixed_summary <- apply(X = fixed_samps,MARGIN = 1, summary)
    colnames(fixed_summary) <- names(object$fixed_samp_indexes)
    fixed_sd <- apply(X = fixed_samps, MARGIN = 1, sd)
    fixed_summary <- rbind(fixed_summary, fixed_sd)
    rownames(fixed_summary)[nrow(fixed_summary)] <- "sd"
    cat("Here are some moments and quantiles for the fixed effects: \n\n")
    print(t(fixed_summary[c(2:5, 7), , drop = FALSE]))
  }
}

#' To predict the GP component in the fitted model, at the locations specified in `newdata`. 
#' @param object The fitted object from the function `model_fit`.
#' @param newdata The dataset that contains the locations to be predicted for the specified GP. Its column names must include `variable`.
#' @param variable The name of the variable to be predicted, should be in the `newdata`.
#' @param degree The degree of derivative that the user specifies for inference. Only applicable for a GP in the `IWP` type.
#' @param include.intercept A logical variable specifying whether the intercept should be accounted when doing the prediction. The default is TRUE. For Coxph model, this 
#' variable will be forced to FALSE.
#' @param only.samples A logical variable indicating whether only the posterior samples are required. The default is FALSE, and the summary of posterior samples will be reported.
#' @export
predict.FitResult <- function(object, newdata = NULL, variable, degree = 0, include.intercept = TRUE, only.samples = FALSE) {
  if(object$family == "Coxph" || object$family == "coxph"| object$family == "cc" | object$family == "casecrossover" | object$family == "CaseCrossover"){
    include.intercept = FALSE ## No intercept for coxph model
  }
  # samps <- aghq::sample_marginal(object$mod, M = 3000)
  samps <- object$samps
  for (instance in object$instances) {
    if(sum(names(object$random_samp_indexes) == variable) >= 2){
      stop("There are more than one variables with the name `variable`: please refit the model with different names for them.")
    }else if(sum(names(object$random_samp_indexes) == variable) == 0){
      stop("The specified variable cannot be found in the fitted model, please check the name.")
    }
    global_samps <- samps$samps[object$boundary_samp_indexes[[variable]], , drop = F]
    coefsamps <- samps$samps[object$random_samp_indexes[[variable]], ]
    if (instance@smoothing_var == variable && class(instance) == "IWP") {
      IWP <- instance
      ## Step 2: Initialization
      if (is.null(newdata)) {
        refined_x_final <- IWP@observed_x
      } else {
        refined_x_final <- sort(newdata[[variable]] - IWP@initial_location) # initialize according to `initial_location`
      }
      if(include.intercept){
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["intercept"]], , drop = F]
      } else{
        intercept_samps <- NULL
      }
      ## Step 3: apply `compute_post_fun_IWP` to samps
      f <- compute_post_fun_IWP(
        samps = coefsamps, global_samps = global_samps,
        knots = IWP@knots,
        refined_x = refined_x_final,
        p = IWP@order, ## check this order or p?
        degree = degree,
        intercept_samps = intercept_samps
      )
      f[,1] <- f[,1] + IWP@initial_location
    }
    else if (instance@smoothing_var == variable && class(instance) == "sGP") {
      sGP <- instance
      ## Step 2: Initialization
      if (is.null(newdata)) {
        refined_x_final <- sGP@observed_x
      } else {
        refined_x_final <- sort(newdata[[variable]] - sGP@initial_location) # initialize according to `initial_location`
      }
      if(include.intercept){
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["intercept"]], , drop = F]
      } else{
        intercept_samps <- NULL
      }
      ## Step 3: apply `compute_post_fun_sGP` to samps
      f <- compute_post_fun_sGP(
        samps = coefsamps, global_samps = global_samps,
        k = sGP@k,
        refined_x = refined_x_final,
        a = sGP@a, 
        m = sGP@m,
        region = sGP@region,
        intercept_samps = intercept_samps,
        boundary = sGP@boundary
      )
      f[,1] <- f[,1] + sGP@initial_location
    }
  }
  if(only.samples){
    return(f)
  }
  ## Step 4: summarize the prediction
  fpos <- extract_mean_interval_given_samps(f)
  names(fpos)[1] <- variable
  return(fpos)
}

#' @export
plot.FitResult <- function(object) {
  ### Step 1: predict with newdata = NULL
  for (instance in object$instances) { ## for each variable in model_fit
    if (class(instance) == "IWP") {
      predict_result <- predict(object, variable = as.character(instance@smoothing_var))
      matplot(
        x = predict_result[,1], y = predict_result[, c("mean", "plower", "pupper")], lty = c(1, 2, 2), lwd = c(2, 1, 1),
        col = "black", type = "l",
        ylab = "effect", xlab = as.character(instance@smoothing_var)
      )
    }
    if (class(instance) == "sGP") {
      predict_result <- predict(object, variable = as.character(instance@smoothing_var))
      matplot(
        x = predict_result[,1], y = predict_result[, c("mean", "plower", "pupper")], lty = c(1, 2, 2), lwd = c(2, 1, 1),
        col = "black", type = "l",
        ylab = "effect", xlab = as.character(instance@smoothing_var)
      )
    }
  }

  ### Step 2: then plot the original aghq object
  # plot(object$mod)
}


#' Extract the posterior samples from the fitted model for the target fixed variables.
#' 
#' @param model_fit The result from model_fit().
#' @param variables A vector of names of the target fixed variables to sample.
#' @export
sample_fixed_effect <- function(model_fit, variables){
  samps <- model_fit$samps$samps
  index <- model_fit$fixed_samp_indexes[variables]
  selected_samps <- t(samps[unlist(index), ,drop = FALSE])
  colnames(selected_samps) <- variables
  return(selected_samps)
} 


#' Computing the posterior samples of the function or its derivative using the posterior samples
#' of the basis coefficients for IWP
#'
#' @param samps A matrix that consists of posterior samples for the O-spline basis coefficients. Each column
#' represents a particular sample of coefficients, and each row is associated with one basis function. This can
#' be extracted using `sample_marginal` function from `aghq` package.
#' @param global_samps A matrix that consists of posterior samples for the global basis coefficients. If NULL,
#' assume there will be no global polynomials and the boundary conditions are exactly zero.
#' @param intercept_samps A matrix that consists of posterior samples for the intercept parameter. If NULL, assume
#' the function evaluated at zero is zero.
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param degree The order of the derivative to take, if zero, implies to consider the function itself.
#' @return A data.frame that contains different samples of the function or its derivative, with the first column
#' being the locations of evaluations x = refined_x.
#' @examples
#' knots <- c(0, 0.2, 0.4, 0.6)
#' samps <- matrix(rnorm(n = (3 * 10)), ncol = 10)
#' result <- compute_post_fun_IWP(samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' global_samps <- matrix(rnorm(n = (2 * 10), sd = 0.1), ncol = 10)
#' result <- compute_post_fun_IWP(global_samps = global_samps, samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' @export
compute_post_fun_IWP <- function(samps, global_samps = NULL, knots, refined_x, p, degree = 0, intercept_samps = NULL) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps <- matrix(0, nrow = (p - 1), ncol = ncol(samps))
  }
  if (nrow(global_samps) != (p - 1)) {
    return(message("Error: Incorrect dimension of global_samps. Check whether the choice of p is consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  if (is.null(intercept_samps)) {
    intercept_samps <- matrix(0, nrow = 1, ncol = ncol(samps))
  }
  if (nrow(intercept_samps) != (1)) {
    return(message("Error: Incorrect dimension of intercept_samps."))
  }
  if (ncol(samps) != ncol(intercept_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and the intercept."))
  }

  ### Augment the global_samps to also consider the intercept
  global_samps <- rbind(intercept_samps, global_samps)

  ## Design matrix for the spline basis weights
  B <- dgTMatrix_wrapper(local_poly_helper(knots, refined_x = refined_x, p = (p - degree)))

  if ((p - degree) >= 1) {
    X <- global_poly_helper(refined_x, p = p)
    X <- as.matrix(X[, 1:(p - degree), drop = FALSE])
    for (i in 1:ncol(X)) {
      X[, i] <- (factorial(i + degree - 1) / factorial(i - 1)) * X[, i]
    }
    fitted_samps_deriv <- X %*% global_samps[(1 + degree):(p), , drop = FALSE] + B %*% samps
  } else {
    fitted_samps_deriv <- B %*% samps
  }
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}



#' Computing the posterior samples of the function using the posterior samples
#' of the basis coefficients for sGP
#'
#' @param samps A matrix that consists of posterior samples for the O-spline basis coefficients. Each column
#' represents a particular sample of coefficients, and each row is associated with one basis function. This can
#' be extracted using `sample_marginal` function from `aghq` package.
#' @param global_samps A matrix that consists of posterior samples for the global basis coefficients. If NULL,
#' assume there will be no global polynomials and the boundary conditions are exactly zero.
#' @param k The number of the sB basis.
#' @param region The region to define the sB basis
#' @param refined_x A vector of locations to evaluate the sB basis
#' @param a The frequency of sGP.
#' @param m The number of harmonics to consider
#' @return A data.frame that contains different samples of the function, with the first column
#' being the locations of evaluations x = refined_x.
#' @export
compute_post_fun_sGP <- function(samps, global_samps = NULL, k, refined_x, a, region, boundary = TRUE, m, intercept_samps = NULL) {
  ## Design matrix for the spline basis weights
  B <- dgTMatrix_wrapper(Compute_B_sB_helper(refined_x = refined_x, k = k, a = a, region = region, boundary = boundary, initial_location = NULL, m = m))
  X <- cbind(1,global_poly_helper_sGP(refined_x = refined_x, a = a, m = m))
  if (is.null(intercept_samps)) {
    intercept_samps <- matrix(0, nrow = 1, ncol = ncol(samps))
  }
  if(is.null(global_samps)){
    global_samps <- matrix(0, nrow = (2*m), ncol = ncol(samps))
  }
  global_samps <- rbind(intercept_samps, global_samps)
  f_samps <- X %*% global_samps + B %*% samps

  result <- cbind(x = refined_x, data.frame(as.matrix(f_samps)))
  result
}


#' Construct posterior inference given samples
#'
#' @param samps Posterior samples of f or its derivative, with the first column being evaluation
#' points x. This can be yielded by `compute_post_fun_IWP` function.
#' @param level The level to compute the pointwise interval.
#' @return A dataframe with a column for evaluation locations x, and posterior mean and pointwise
#' intervals at that set of locations.
#' @export
extract_mean_interval_given_samps <- function(samps, level = 0.95) {
  x <- samps[, 1]
  samples <- samps[, -1]
  result <- data.frame(x = x)
  alpha <- 1 - level
  result$plower <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (alpha / 2)))
  result$pupper <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (level + (alpha / 2))))
  result$mean <- as.numeric(apply(samples, MARGIN = 1, mean))
  result
}




#' Obtain the posterior density of a variance parameter in the fitted model
#' 
#' @param object The fitted object from the function `model_fit`.
#' @param component The component of the variance parameter that you want to show. By default this is `NULL`, indicating the family.sd is of interest.
#' @param h For PSD, the unit of predictive step to consider, by default is set to `NULL`, indicating the result is using the same `h` as in the model fitting.
#' @param theta_logprior The log prior function used on the selected variance parameter. By default is `NULL`, and the default Exponential prior will be used.
#' @param MCMC_samps_only For model fitted with MCMC, whether only the posterior samples are needed.
#' @export
var_density <- function(object, component = NULL, h = NULL, theta_logprior = NULL, MCMC_samps_only = FALSE){
  postsigma <- NULL
  
  if(is.null(theta_logprior)){
    theta_logprior <- function(theta,prior_alpha, prior_u) {
      lambda <- -log(prior_alpha)/prior_u
      log(lambda/2) - lambda * exp(-theta/2) - theta/2
    }
  }
  priorfunc <- function(x,prior_alpha, prior_u) exp(theta_logprior(x,prior_alpha, prior_u))
  priorfuncsigma <- function(x,prior_alpha, prior_u) (2/x) * exp(theta_logprior(-2*log(x), prior_alpha, prior_u))
  
  if(any(class(object$mod) == "aghq")){
    if(is.null(component)){
      if(object$family != "Gaussian"){
        stop("There is no family SD in the fitted model. Please indicate which component of the var-parameter that you want to show in `component`.")
      }
      theta_marg <- object$mod$marginals[[length(object$instances) + 1]]
      if(nrow(theta_marg) <= 2){
        stop("The number of quadrature points is too small, please use aghq_k >= 3.")
      }
      logpostsigma <- compute_pdf_and_cdf(theta_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
      postsigma <- data.frame(SD = logpostsigma$transparam, 
                                  post = logpostsigma$pdf_transparam,
                                  prior = priorfuncsigma(logpostsigma$transparam, prior_alpha = object$control.family$sd.prior$param$alpha, prior_u = object$control.family$sd.prior$param$u))
    }
    else{
      for (i in 1:length(object$instances)) {
        instance <- object$instances[[i]]
        if(instance@smoothing_var == component){
          theta_marg <- object$mod$marginals[[i]]
          if(nrow(theta_marg) <= 2){
            stop("The number of quadrature points is too small, please use aghq_k >= 3.")
          }
          logpostsigma <- compute_pdf_and_cdf(theta_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
          postsigma <- data.frame(SD = logpostsigma$transparam, 
                                  post = logpostsigma$pdf_transparam,
                                  prior = priorfuncsigma(logpostsigma$transparam, prior_alpha = object$instances[[i]]@sd.prior$param$alpha, prior_u = object$instances[[i]]@sd.prior$param$u))
          
          if(is.null(h)){
            if(!is.null(instance@sd.prior$h)){
              h <- instance@sd.prior$h
            }
          }
          if(!is.null(h)){
            if(class(instance) == "IWP"){
              p <- instance@order
              correction <- sqrt((h^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
            }
            else if(class(instance) == "sGP"){
              correction <- 0
              for (j in 1:instance@m) {
                correction <- correction + compute_d_step_sGPsd(d = h, a = (j*instance@a))
              }
            }
            else{
              stop("PSD is currently on defined on IWP and sGP, please specify h = NULL for other type of random effect")
            }
            postsigmaPSD <- data.frame(PSD = postsigma$SD * correction, 
                                       post.PSD = postsigma$post / correction,
                                       prior.PSD = postsigma$prior / correction)
            
            postsigma <- cbind(postsigma, postsigmaPSD)
          }
          
        }
      }
      
    }
  }
  
  else if(any(class(object$mod) == "stanfit")){
    sigmaPSD_marg_samps <- NULL
    if(is.null(component)){
      if(object$family != "Gaussian"){
        stop("There is no family SD in the fitted model. Please indicate which component of the var-parameter that you want to show in `component`.")
      }
      theta_marg_samps <- object$samps$thet[[length(object$instances) + 1]]
      sigma_marg_samps <- exp(-0.5*theta_marg_samps)
      sigma_marg_density <- density(sigma_marg_samps)
      postsigma <- data.frame(SD = sigma_marg_density$x, 
                              post = sigma_marg_density$y,
                              prior = priorfuncsigma(sigma_marg_density$x, prior_alpha = object$control.family$sd.prior$param$alpha, prior_u = object$control.family$sd.prior$param$u))
    }
    else{
      for (i in 1:length(object$instances)) {
        instance <- object$instances[[i]]
        if(instance@smoothing_var == component){
          theta_marg_samps <- object$samps$thet[[i]]
          sigma_marg_samps <- exp(-0.5*theta_marg_samps)
          sigma_marg_density <- density(sigma_marg_samps)
          postsigma <- data.frame(SD = sigma_marg_density$x, 
                                  post = sigma_marg_density$y,
                                  prior = priorfuncsigma(sigma_marg_density$x, prior_alpha = object$instances[[i]]@sd.prior$param$alpha, prior_u = object$instances[[i]]@sd.prior$param$u))
          if(is.null(h)){
            if(!is.null(instance@sd.prior$h)){
              h <- instance@sd.prior$h
            }
          }
          if(!is.null(h)){
            if(class(instance) == "IWP"){
              p <- instance@order
              correction <- sqrt((h^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
            }
            else if(class(instance) == "sGP"){
              correction <- 0
              for (j in 1:instance@m) {
                correction <- correction + compute_d_step_sGPsd(d = h, a = (j*instance@a))
              }
            }
            else{
              stop("PSD is currently on defined on IWP and sGP, please specify h = NULL for other type of random effect")
            }
            sigmaPSD_marg_samps <- sigma_marg_samps * correction
            postsigmaPSD <- data.frame(PSD = postsigma$SD * correction, 
                                       post.PSD = postsigma$post / correction,
                                       prior.PSD = postsigma$prior / correction)
            
            postsigma <- cbind(postsigma, postsigmaPSD)
          }
          
        }
      }
    }
    if(MCMC_samps_only == TRUE){
      return(list(sigmaPSD_marg_samps = sigmaPSD_marg_samps, sigma_marg_samps = sigma_marg_samps))
    }
  }
  
  else{
    stop("The function `var_density` currently only supports model fittd with `method = aghq`.")
  }
  postsigma <- postsigma[order(postsigma$SD),]
  return(postsigma)
}


#' Obtain the posterior and prior density of all the parameters in the fitted model
#' 
#' @param object The fitted object from the function `model_fit`.
#' @export
para_density <- function(object){
  result_list <- list()
  for (fixed_name in names(object$fixed_samp_indexes)) {
    samps <- sample_fixed_effect(object, variables = fixed_name)
    fixed_density <- density(samps)
    result_list[[fixed_name]] <- data.frame(effect = fixed_density$x, post = fixed_density$y)
  }
  
  for (random_name in names(object$random_samp_indexes)) {
    result_list[[random_name]] <- var_density(object = object, component = random_name)
  }
  
  if(object$family == "Gaussian"){
    result_list[["family_var"]] <- var_density(object = object)
  }
  
  return(result_list)
}


#' Obtain the posterior summary table for all the parameters in the fitted model
#' @param quantiles The specified quantile to display the posterior summary, default is c(0.025, 0.975).
#' @param digits The significant digits to be kept in the result, default is 3.
#' @export
post_table <- function(object, quantiles = c(0.025, 0.975), digits = 3){
  all_density <- para_density(object)
  all_cdf <- list()
  compute_cdf <- function(x, y){
    new_data_frame <- list(x = x)
    new_data_frame$cdf <- cumsum(y * c(diff(x), 0))
    new_data_frame
  }
  result_table <- c("name", "median", paste0("q", unlist(strsplit(toString(quantiles), split = ", "))), "prior", "prior:P1", "prior:P2")
  for (name in names(object$fixed_samp_indexes)) {
    all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$effect, y = all_density[[name]]$post)
    to_add <- c(name, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= q))])
    }
    to_add <- c(to_add, "Normal", object$control.fixed[[name]]$mean, 1/object$control.fixed[[name]]$prec)
    result_table <- rbind(result_table, to_add)
  } 
  
  for (name in names(object$random_samp_indexes)) {
    if("PSD" %in% names(all_density[[name]])){
      all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$PSD, y = all_density[[name]]$post.PSD)
      print_name <- paste0(name, " (PSD)")
    }
    else{
      all_cdf[[name]] <- compute_cdf(x = all_density[[name]]$SD, y = all_density[[name]]$post)
      print_name <- paste0(name, " (SD)")
    }
    to_add <- c(print_name, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[[name]]$x[max(which(all_cdf[[name]]$cdf <= q))])
    }
    for (i in 1:length(object$instances)) {
      if(name == object$instances[[i]]@smoothing_var){
        param <- object$instances[[i]]@sd.prior$param
      }
    }
    object$instances
    to_add <- c(to_add, "Exponential", param$u, param$alpha)
    result_table <- rbind(result_table, to_add)
  }
  
  if("family_var" %in% names(all_density)){
    all_cdf[["family_var"]] <- compute_cdf(x = all_density[["family_var"]]$SD, y = all_density[["family_var"]]$post)
    to_add <- c("family_var", all_cdf[["family_var"]]$x[max(which(all_cdf[["family_var"]]$cdf <= 0.5))])
    for (q in quantiles) {
      to_add <- c(to_add, all_cdf[["family_var"]]$x[max(which(all_cdf[["family_var"]]$cdf <= q))])
    }
    to_add <- c(to_add, "Exponential", object$control.family$sd.prior$param$u, object$control.family$sd.prior$param$alpha)
    result_table <- rbind(result_table, to_add)
  }
  result_table <- data.frame(unname(result_table))
  result_table_names <- unlist(unname(result_table)[1,])
  result_table <- (type.convert(result_table[-1,], as.is = TRUE))
  result_table <- data.frame(lapply(result_table, function(y) if(is.numeric(y)) round(y, digits) else y)) 
  colnames(result_table) <- result_table_names
  result_table
}


