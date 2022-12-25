#' Constructing and evaluating the local O-spline basis (design matrix)
#'
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of refined_x, the ncol should equal to number of knots minus 1, and nrow
#' should equal to the number of elements in refined_x.
#' @examples
#' local_poly(knots = c(0,0.2,0.4,0.6,0.8), refined_x = seq(0,0.8, by = 0.1), p = 2)
#' @export
local_poly <- function(knots, refined_x, p = 2){
  dif <- diff(knots)
  nn <- length(refined_x)
  n <- length(knots)
  D <- matrix(0, nrow = nn, ncol = n-1)
  for (j in 1:nn) {
    for (i in 1:(n-1)) {
      if(refined_x[j] <= knots[i]){
        D[j, i] <- 0
      }
      else if (refined_x[j] <= knots[i+1] & refined_x[j] >= knots[i]){
        D[j, i] <- (1/factorial(p)) * (refined_x[j] - knots[i])^p
      }
      else{
        k <- 1:p
        D[j,i] <- sum((dif[i]^k)*((refined_x[j] - knots[i+1])^(p-k))/(factorial(k)*factorial(p-k)))
      }
      }
  }
  D
}

#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix)
#'
#' @param x A vector of locations to evaluate the global polynomials
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to p, and nrow
#' should equal to the number of elements in x
#' @examples
#' global_poly(x = c(0,0.2,0.4,0.6,0.8), p = 2)
#' @export
global_poly <- function(x, p = 2){
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i-1))
  }
  result
}

#' Constructing the precision matrix given the knot sequence
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @examples
#' compute_weights_precision(x = c(0,0.2,0.4,0.6,0.8))
#' @export
compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}

#' Computing the posterior samples of the function or its derivative using the posterior samples
#' of the basis coefficients
#'
#' @param samps A matrix that consists of posterior samples for the O-spline basis coefficients. Each column
#' represents a particular sample of coefficients, and each row is associated with one basis function. This can
#' be extracted using `sample_marginal` function from `aghq` package.
#' @param global_samps A matrix that consists of posterior samples for the global basis coefficients. If NULL,
#' assume there will be no global polynomials and the boundary conditions are exactly zero.
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @param degree The order of the derivative to take, if zero, implies to consider the function itself.
#' @return A data.frame that contains different samples of the function or its derivative, with the first column
#' being the locations of evaluations x = refined_x.
#' @examples
#' knots <- c(0,0.2,0.4,0.6)
#' samps <- matrix(rnorm(n = (3*10)), ncol = 10)
#' result <- compute_post_fun(samps = samps, knots = knots, refined_x = seq(0,1, by = 0.1), p = 2)
#' plot(result[,2]~result$x, type = 'l', ylim = c(-0.3,0.3))
#' for (i in 1:9) {
#'   lines(result[,(i+1)]~result$x, lty = 'dashed', ylim = c(-0.1,0.1))
#' }
#' global_samps <- matrix(rnorm(n = (2*10), sd = 0.1), ncol = 10)
#' result <- compute_post_fun(global_samps = global_samps, samps = samps, knots = knots, refined_x = seq(0,1, by = 0.1), p = 2)
#' plot(result[,2]~result$x, type = 'l', ylim = c(-0.3,0.3))
#' for (i in 1:9) {
#'   lines(result[,(i+1)]~result$x, lty = 'dashed', ylim = c(-0.1,0.1))
#' }
#' @export
compute_post_fun <- function(samps, global_samps = NULL, knots, refined_x, p, degree = 0){
  if(p <= degree){return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))}
  if(is.null(global_samps)){global_samps = matrix(0, nrow = p, ncol = ncol(samps))}
  if(nrow(global_samps)!= p | nrow(samps)!= (length(knots)-1)){
    return(message("Error: Incorrect dimension of global_samps or samps. Check whether the choice of p or the choice of knots are consistent with the fitted model."))
  }
  if(ncol(samps)!=ncol(global_samps)){
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  X = global_poly(refined_x, p = p)
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  ## Design matrix for the spline basis weights
  B = as(local_poly(knots, refined_x = refined_x, p = (p-degree)),"dgTMatrix")
  fitted_samps_deriv <- X %*% global_samps[(1 + degree):p, ] + B %*% samps
  result <- cbind(x = refined_x, data.frame(as.matrix(fitted_samps_deriv)))
  result
}


#' Construct posterior inference given samples
#'
#' @param samps Posterior samples of f or its derivative, with the first column being evaluation
#' points x. This can be yielded by `compute_post_fun` function.
#' @param level The level to compute the pointwise interval.
#' @return A dataframe with a column for evaluation locations x, and posterior mean and pointwise
#' intervals at that set of locations.
#' @export
extract_mean_interval_given_samps <- function(samps, level = 0.95){
  x <- samps[,1]
  samples <- samps[,-1]
  result <- data.frame(x = x)
  alpha <- 1-level
  result$plower <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (alpha/2)))
  result$pupper <- as.numeric(apply(samples, MARGIN = 1, quantile, p = (level + (alpha/2))))
  result$mean <- as.numeric(apply(samples, MARGIN = 1, mean))
  result
}



#' Construct prior based on d-step prediction SD.
#'
#' @param prior A list that contains a and u. This specifies the target prior on the d-step SD \equ{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = a}.
#' @param d A numeric value for the prediction step.
#' @param p An integer for the order of IWP.
#' @return A list that contains a and u. The prior for the smoothness parameter \equ{\sigma} such that \eqn{P(\sigma > u) = a}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion <- function(d, prior, p){
  Cp <- (d^((2*p) - 1))/(((2*p)-1)*(factorial(p-1)^2))
  prior_q <- list(a = prior$a, u = (prior$u * (1 /sqrt(Cp))))
  prior_q
}






#' Roxygen commands
#'
#' @useDynLib OSplines
#'
dummy <- function(){
  return(NULL)
}
