#' @export
summary.FitResult <- function(object) {
  if(any(class(object$mod) == "aghq")){
    aghq_summary <- summary(object$mod)
    aghq_output <- capture.output(aghq_summary)
    cur_index <- grep("Here are some moments and quantiles for the transformed parameter:", aghq_output)
    cat(paste(aghq_output[1:(cur_index - 1)], collapse = "\n"))
    cat("\nHere are some moments and quantiles for the log precision: \n")
    cur_index <- grep("median    97.5%", aghq_output)
    cat(paste(aghq_output[cur_index], collapse = "\n"))
    
    summary_table <- as.matrix(aghq_summary$summarytable)
    theta_names <- c()
    for (instance in object$instances) {
      theta_names <- c(theta_names, paste("theta(", toString(instance@smoothing_var), ")", sep = ""))
    }
    
    row.names(summary_table) <- theta_names
    print(summary_table)
    cat("\n")
  }
  # samps <- aghq::sample_marginal(object$mod, M = 3000)
  samps <- object$samps
  fixed_samps <- samps$samps[unlist(object$fixed_samp_indexes), , drop = F]
  fixed_summary <- apply(X = fixed_samps,MARGIN = 1, summary)
  colnames(fixed_summary) <- names(object$fixed_samp_indexes)
  fixed_sd <- apply(X = fixed_samps, MARGIN = 1, sd)
  fixed_summary <- rbind(fixed_summary, fixed_sd)
  rownames(fixed_summary)[nrow(fixed_summary)] <- "sd"

  cat("Here are some moments and quantiles for the fixed effects: \n\n")
  print(t(fixed_summary[c(2:5, 7), ]))
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
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["Intercept"]], , drop = F]
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
        intercept_samps <- samps$samps[object$fixed_samp_indexes[["Intercept"]], , drop = F]
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
sample_fixed_effect <- function(model_fit, variables, M){
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
