# Call with ::
library(tidyverse)
library(Matrix)
library(TMB)
library(aghq)

#' @title
#' Model fitting with random effects/fixed effects
#'
#' @description
#' Fitting a hierarchical model based on the provided formula, data and parameters such as type of method and family of response.
#' Returning the S4 objects for the random effects, concatenated design matrix for the intercepts and fixed effects, fitted model,
#' indexes to partition the posterior samples.
#'
#' @param formula A formula that contains one response variable, and covariates with either random or fixed effect.
#' @param data A dataframe that contains the response variable and other covariates mentioned in the formula.
#' @param method The inference method used in the model. By default, the method is set to be "aghq".
#' @param family The family of response used in the model. By default, the family is set to be "Gaussian".
#' @param control.family Parameters used to specify the priors for the family parameters, such as the standard deviation parameter of Gaussian family.
#' @param control.fixed Parameters used to specify the priors for the fixed effects.
#' @return A list that contains following items: the S4 objects for the random effects (instances), concatenated design matrix for
#' the fixed effects (design_mat_fixed), fitted aghq (mod) and indexes to partition the posterior samples
#' (boundary_samp_indexes, random_samp_indexes and fixed_samp_indexes).
#' @examples
#' library(OSplines)
#' library(tidyverse)
#'
#' data <- INLA::Munich %>% select(rent, floor.size, year, location)
#' data$score <- rnorm(n = nrow(data))
#' head(data, n = 5)
#'
#' ### A model with two IWP and two Fixed effects:
#' ## Assume f(floor.size) is second order IWP
#' ## Assume f(year) is third order IWP
#'
#' fit_result <- model_fit(
#'   rent ~ location + f(
#'     smoothing_var = floor.size,
#'     model = "IWP",
#'     order = 2
#'   )
#'   + score + f(
#'       smoothing_var = year,
#'       model = "IWP",
#'       order = 3, k = 10, # should add a checker for k >= 3
#'       sd.prior = list(prior = "exp", para = list(u = 1, alpha = 0.5)),
#'       boundary.prior = list(prec = 0.01)
#'     ),
#'   data = data, method = "aghq", family = "Gaussian",
#'   control.family = list(sd_prior = list(prior = "exp", para = list(u = 1, alpha = 0.5))),
#'   control.fixed = list(intercept = list(prec = 0.01), location = list(prec = 0.01), score = list(prec = 0.01))
#' )
#'
#' # Check the contents of the returned fit result
#' names(fit_result)
#' IWP1 <- fit_result$instances[[1]]
#' IWP2 <- fit_result$instances[[2]]
#' mod <- fit_result$mod
#' @export
model_fit <- function(formula, data, method = "aghq", family = "Gaussian", control.family, control.fixed) {
  # parse the input formula
  parse_result <- parse_formula(formula)
  response_var <- parse_result$response
  rand_effects <- parse_result$rand_effects
  fixed_effects <- parse_result$fixed_effects

  instances <- list()
  design_mat_fixed <- list()

  # For random effects
  for (rand_effect in rand_effects) {
    smoothing_var <- rand_effect$smoothing_var
    model_class <- rand_effect$model
    sd.prior <- eval(rand_effect$sd.prior)
    if (is.null(sd.prior)) {
      sd.prior <- list(prior = "exp", para = list(u = 1, alpha = 0.5))
    }
    if (sd.prior$prior != "exp") {
      stop("Error: For each random effect, sd.prior currently only supports 'exp' (exponential) as prior.")
    }
    cat(model_class)
    if(model_class == "IWP"){
      order <- eval(rand_effect$order)
      knots <- eval(rand_effect$knots)
      k <- eval(rand_effect$k)
      initial_location <- eval(rand_effect$initial_location)
      if (!(is.null(k)) && k < 3) {
        stop("Error: Parameter <k> in the random effect part should be >= 3.")
      }
      if(order < 1){
        stop("Error: Parameter <order> in the random effect part should be >= 1.")
      }
      boundary.prior <- eval(rand_effect$boundary.prior)
      # If the user does not specify initial_location, compute initial_location with
      # the min of data[[smoothing_var]]
      if (is.null(initial_location)){
        initial_location = min(data[[smoothing_var]])
      }
      # If the user does not specify knots, compute knots with
      # the parameter k
      if (is.null(knots)) {
        initialized_smoothing_var <- data[[smoothing_var]] - initial_location
        default_k <- 5
        if (is.null(k)) {
          knots <- unique(sort(seq(from = min(initialized_smoothing_var), to = max(initialized_smoothing_var), length.out = default_k))) # should be length.out
        } else {
          knots <- unique(sort(seq(from = min(initialized_smoothing_var), to = max(initialized_smoothing_var), length.out = k)))
        }
      }
      # refined_x <- seq(from = min(initialized_smoothing_var), to = max(initialized_smoothing_var), by = 1) # this is not correct
      observed_x <- sort(initialized_smoothing_var) # initialized_smoothing_var: initialized observed covariate values
      if (is.null(boundary.prior)) {
        boundary.prior <- list(prec = 0.01)
      }
      instance <- new(model_class,
                      response_var = response_var,
                      smoothing_var = smoothing_var, order = order,
                      knots = knots, observed_x = observed_x, sd.prior = sd.prior, boundary.prior = boundary.prior, data = data
      )
      # Case for IWP
      instance@initial_location <- initial_location
      instance@X <- global_poly(instance)[, -1, drop = FALSE]
      instance@B <- local_poly(instance)
      instance@P <- compute_weights_precision(instance)
      instances[[length(instances) + 1]] <- instance
    }
    else if(model_class == "IID"){
      instance <- new(model_class,
                      response_var = response_var,
                      smoothing_var = smoothing_var, sd.prior = sd.prior, data = data
      )
      # Case for IID
      instance@B <- compute_B(instance)
      instance@P <- compute_P(instance)
      instances[[length(instances) + 1]] <- instance
    }
  }

  # For the intercept
  Xf0 <- matrix(1, nrow = nrow(data), ncol = 1)
  design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf0


  fixed_effects_names <- c("Intercept")
  # For fixed effects
  for (fixed_effect in fixed_effects) {
    Xf <- matrix(data[[fixed_effect]], nrow = nrow(data), ncol = 1)
    design_mat_fixed[[length(design_mat_fixed) + 1]] <- Xf
    fixed_effects_names <- c(fixed_effects_names, fixed_effect)
  }

  if (missing(control.family)) {
    control.family <- list(sd_prior = list(prior = "exp", para = list(u = 1, alpha = 0.5)))
  }

  if (control.family$sd_prior$prior != "exp") {
    stop("Error: Currently, control.family only supports 'exp' (exponential) as prior.")
  }

  if (missing(control.fixed)) {
    control.fixed <- list(intercept = list(prec = 0.01))
    for (fixed_effect in fixed_effects) {
      control.fixed$fixed_effect <- list(prec = 0.01)
    }
  }

  result_by_method <- get_result_by_method(instances, design_mat_fixed, family, control.family, control.fixed, fixed_effects)
  mod <- result_by_method$mod
  w_count <- result_by_method$w_count

  global_samp_indexes <- list()
  coef_samp_indexes <- list()
  fixed_samp_indexes <- list()
  rand_effects_names <- c()
  global_effects_names <- c()
  sum_col_ins <- 0
  for (instance in instances) {
    sum_col_ins <- sum_col_ins + ncol(instance@B)
    rand_effects_names <- c(rand_effects_names, instance@smoothing_var)
    if(class(instance) == "IWP"){
      global_effects_names <- c(global_effects_names, instance@smoothing_var)
    }
  }

  cur_start <- sum_col_ins + 1
  cur_end <- sum_col_ins
  cur_coef_start <- 1
  cur_coef_end <- 0
  for (instance in instances) {
    if(class(instance) == "IWP"){
      cur_end <- cur_end + ncol(instance@X)
      if(instance@order == 1){
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- numeric()
      }
      else if(instance@order > 1){
        global_samp_indexes[[length(global_samp_indexes) + 1]] <- (cur_start:cur_end)
      }
    }
    cur_coef_end <- cur_coef_end + ncol(instance@B)
    coef_samp_indexes[[length(coef_samp_indexes) + 1]] <- (cur_coef_start:cur_coef_end)
    cur_start <- cur_end + 1
    cur_coef_start <- cur_coef_end + 1
  }
  names(global_samp_indexes) <- global_effects_names
  names(coef_samp_indexes) <- rand_effects_names

  for (fixed_samp_index in ((cur_end + 1):w_count)) {
    fixed_samp_indexes[[length(fixed_samp_indexes) + 1]] <- fixed_samp_index
  }

  names(fixed_samp_indexes) <- fixed_effects_names

  return(list(
    instances = instances, design_mat_fixed = design_mat_fixed, mod = mod,
    boundary_samp_indexes = global_samp_indexes,
    random_samp_indexes = coef_samp_indexes,
    fixed_samp_indexes = fixed_samp_indexes
  ))
}

parse_formula <- function(formula) {
  components <- as.list(attributes(terms(formula))$ variables)
  fixed_effects <- list()
  rand_effects <- list()
  # Index starts as 3 since index 1 represents "list" and
  # index 2 represents the response variable
  for (i in 3:length(components)) {
    if (startsWith(toString(components[[i]]), "f,")) {
      rand_effects[[length(rand_effects) + 1]] <- components[[i]]
    } else {
      fixed_effects[[length(fixed_effects) + 1]] <- components[[i]]
    }
  }
  return(list(response = components[[2]], fixed_effects = fixed_effects, rand_effects = rand_effects))
}

# Create a class for IWP using S4
setClass("IWP", slots = list(
  response_var = "name", smoothing_var = "name", order = "numeric",
  knots = "numeric", observed_x = "numeric", sd.prior = "list",
  boundary.prior = "list", data = "data.frame", X = "matrix",
  B = "matrix", P = "matrix", initial_location = "numeric"
))

# Create a class for IID using S4
setClass("IID", slots = list(
  response_var = "name", smoothing_var = "name", sd.prior = "list",
  data = "data.frame", B = "matrix", P = "matrix"
))

setGeneric("compute_B", function(object) {
  standardGeneric("compute_B")
})
setMethod("compute_B", signature = "IID", function(object) {
  smoothing_var <- object@smoothing_var
  x <- as.factor((object@data)[[smoothing_var]])
  B <- model.matrix(~-1+x)
  B
})

setGeneric("compute_P", function(object) {
  standardGeneric("compute_P")
})
setMethod("compute_P", signature = "IID", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]]
  num_factor <- length(unique(x))
  diag(nrow = num_factor, ncol = num_factor)
})

setGeneric("local_poly", function(object) {
  standardGeneric("local_poly")
})
setMethod("local_poly", signature = "IWP", function(object) {
  knots <- object@knots
  initial_location <- object@initial_location
  smoothing_var <- object@smoothing_var
  refined_x <- (object@data)[[smoothing_var]] - initial_location
  p <- object@order
  # TODO: refactor this part
  if (min(knots) >= 0){ 
    # The following part only works with all-positive knots
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= knots[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x[j] - knots[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - knots[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
  }
  else if(max(knots) <= 0){
    # Handle the negative part only
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] & refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] - knots_neg[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - knots_neg[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
  }
  else{
    # Handle the negative part
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] & refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] - knots_neg[i])^p
        } else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - knots_neg[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
    
    # Handle the positive part
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        } else if (refined_x_pos[j] <= knots_pos[i + 1] & refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1 / factorial(p)) * (refined_x_pos[j] - knots_pos[i])^p
        } else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] - knots_pos[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D
})

setGeneric("global_poly", function(object) {
  standardGeneric("global_poly")
})
setMethod("global_poly", signature = "IWP", function(object) {
  smoothing_var <- object@smoothing_var
  x <- (object@data)[[smoothing_var]] - object@initial_location
  p <- object@order
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
})


#' Constructing the precision matrix given the knot sequence
#'
#' @param x A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @return A precision matrix of the corresponding basis function, should be diagonal matrix with
#' size (k-1) by (k-1).
#' @export
setGeneric("compute_weights_precision", function(object) {
  standardGeneric("compute_weights_precision")
})
setMethod("compute_weights_precision", signature = "IWP", function(object) {
  knots <- object@knots
  if(min(knots) >= 0){
    as(diag(diff(knots)), "matrix")
  }
  else if(max(knots) < 0){
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    as(diag(diff(knots_neg)), "matrix")
  }
  else{
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    d1 <- diff(knots_neg)
    d2 <- diff(knots_pos)
    Precweights1 <- diag(d1)
    Precweights2 <- diag(d2)
    as(Matrix::bdiag(Precweights1, Precweights2), "matrix")
    }
})


get_result_by_method <- function(instances, design_mat_fixed, family, control.family, control.fixed, fixed_effects, aghq_k = 4, size = NULL) {
  # Family types: Gaussian - 0, Poisson - 1, Binomial - 2
  if (family == "Gaussian"){
    family_type = 0
  }
  else if (family == "Poisson"){
    family_type = 1
  }
  else if (family == "Binomial"){
    family_type = 2
  }
  
  # Containers for random effects
  X <- list()
  B <- list()
  P <- list()
  logPdet <- list()
  u <- list()
  alpha <- list()
  betaprec <- list()

  # Containers for fixed effects
  beta_fixed_prec <- list()
  Xf <- list()

  w_count <- 0
  # Need a theta for the Gaussian variance, so
  # theta_count starts at 1 if Gaussian
  theta_count <- 0 + (family_type == 0)

  for (instance in instances) {
    # For each random effects
    if(class(instance) == "IWP"){
      X[[length(X) + 1]] <- dgTMatrix_wrapper(instance@X)
      betaprec[[length(betaprec) + 1]] <- instance@boundary.prior$prec
    }
    B[[length(B) + 1]] <- dgTMatrix_wrapper(instance@B)
    P[[length(P) + 1]] <- dgTMatrix_wrapper(instance@P)
    logPdet[[length(logPdet) + 1]] <- as.numeric(determinant(instance@P, logarithm = TRUE)$modulus)
    u[[length(u) + 1]] <- instance@sd.prior$para$u
    alpha[[length(alpha) + 1]] <- instance@sd.prior$para$alpha
    if(class(instance) == "IWP"){
      w_count <- w_count + ncol(instance@X)
    }
    w_count <- w_count + ncol(instance@B)
    theta_count <- theta_count + 1
  }

  # For the variance of the Gaussian family
  # From control.family, if applicable
  if (family_type == 0){
    u[[length(u) + 1]] <- control.family$sd_prior$para$u
    alpha[[length(alpha) + 1]] <- control.family$sd_prior$para$alpha
  }
  for (i in 1:length(design_mat_fixed)) {
    # For each fixed effects
    if (i == 1) {
      beta_fixed_prec[[i]] <- control.fixed$intercept$prec
    } else {
      beta_fixed_prec[[i]] <- control.fixed[[fixed_effects[[i - 1]]]]$prec
    }
    Xf[[length(Xf) + 1]] <- dgTMatrix_wrapper(design_mat_fixed[[i]])
    w_count <- w_count + ncol(design_mat_fixed[[i]])
  }

  tmbdat <- list(
    # For Random effects
    X = X,
    B = B,
    P = P,
    logPdet = logPdet,
    u = u,
    alpha = alpha,
    betaprec = betaprec,

    # For Fixed Effects:
    beta_fixed_prec = beta_fixed_prec,
    Xf = Xf,

    # Response
    y = (instances[[1]]@data)[[instances[[1]]@response_var]],

    # Family type
    family_type = family_type
  )
  
  # If Family == "Binomial", check whether size is defined in user's input
  if(family_type == 2 & is.null(size)){
    tmbdat$size <- numeric(length = length(tmbdat$y)) + 1 # A vector of 1s being default
  }

  tmbparams <- list(
    W = c(rep(0, w_count)), # recall W is everything in the model (RE or FE)
    theta = c(rep(0, theta_count))
  )

  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "OSplines",
    silent = TRUE
  )

  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
  mod <- aghq::marginal_laplace_tmb(ff, aghq_k, c(rep(0, theta_count))) # The default value of aghq_k is 4

  return(list(mod = mod, w_count = w_count))
}


#' Constructing and evaluating the local O-spline basis (design matrix)
#'
#' @param knots A vector of knots used to construct the O-spline basis, first knot should be viewed as "0",
#' the reference starting location. These k knots will define (k-1) basis function in total.
#' @param refined_x A vector of locations to evaluate the O-spline basis
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j component being the value of jth basis function
#' value at ith element of refined_x, the ncol should equal to number of knots minus 1, and nrow
#' should equal to the number of elements in refined_x.
#' @examples
#' local_poly(knots = c(0, 0.2, 0.4, 0.6, 0.8), refined_x = seq(0, 0.8, by = 0.1), p = 2)
#' @export
local_poly_helper <- function(knots, refined_x, p = 2) {
  # TODO: refactor this part
  if (min(knots) >= 0){ 
    # The following part only works with all-positive knots
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        } else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= knots[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x[j] - knots[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - knots[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
  }
  else if(max(knots) <= 0){
    # Handle the negative part only
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] & refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] - knots_neg[i])^p
        } else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - knots_neg[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
  }
  else{
    # Handle the negative part
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        } else if (refined_x_neg[j] <= knots_neg[i + 1] & refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1 / factorial(p)) * (refined_x_neg[j] - knots_neg[i])^p
        } else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - knots_neg[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
    
    # Handle the positive part
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        } else if (refined_x_pos[j] <= knots_pos[i + 1] & refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1 / factorial(p)) * (refined_x_pos[j] - knots_pos[i])^p
        } else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] - knots_pos[i + 1])^(p - k)) / (factorial(k) * factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D # Local poly design matrix
}


#' Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix)
#'
#' @param x A vector of locations to evaluate the global polynomials
#' @param p An integer value indicates the order of smoothness
#' @return A matrix with i,j componet being the value of jth basis function
#' value at ith element of x, the ncol should equal to p, and nrow
#' should equal to the number of elements in x
#' @examples
#' global_poly(x = c(0, 0.2, 0.4, 0.6, 0.8), p = 2)
#' @export
global_poly_helper <- function(x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
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
#' knots <- c(0, 0.2, 0.4, 0.6)
#' samps <- matrix(rnorm(n = (3 * 10)), ncol = 10)
#' result <- compute_post_fun(samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' global_samps <- matrix(rnorm(n = (2 * 10), sd = 0.1), ncol = 10)
#' result <- compute_post_fun(global_samps = global_samps, samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
#' plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
#' for (i in 1:9) {
#'   lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
#' }
#' @export
compute_post_fun <- function(samps, global_samps = NULL, knots, refined_x, p, degree = 0) {
  if (p <= degree) {
    return(message("Error: The degree of derivative to compute is not defined. Should consider higher order smoothing model or lower order of the derivative degree."))
  }
  if (is.null(global_samps)) {
    global_samps <- matrix(0, nrow = (p-1), ncol = ncol(samps))
  }
  if (nrow(global_samps) != (p-1)) {
    return(message("Error: Incorrect dimension of global_samps. Check whether the choice of p is consistent with the fitted model."))
  }
  if (ncol(samps) != ncol(global_samps)) {
    return(message("Error: The numbers of posterior samples do not match between the O-splines and global polynomials."))
  }
  ## Design matrix for the spline basis weights
  B <- dgTMatrix_wrapper(local_poly_helper(knots, refined_x = refined_x, p = (p - degree)))
  
  if ((p - degree) > 1){
    X <- global_poly_helper(refined_x, p = p)
    X <- as.matrix(X[, 1:(p - degree)])
    for (i in 1:ncol(X)) {
      X[, i] <- (factorial(i + degree - 1) / factorial(i - 1)) * X[, i]
    }
    fitted_samps_deriv <- X[,-1] %*% global_samps[(1 + degree):(p-1), ,drop = FALSE] + B %*% samps
  }
  else {
    fitted_samps_deriv <- B %*% samps
  }
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


#' Construct prior based on d-step prediction SD.
#'
#' @param prior A list that contains a and u. This specifies the target prior on the d-step SD \eqn{\sigma(d)}, such that \eqn{P(\sigma(d) > u) = a}.
#' @param d A numeric value for the prediction step.
#' @param p An integer for the order of IWP.
#' @return A list that contains a and u. The prior for the smoothness parameter \eqn{\sigma} such that \eqn{P(\sigma > u) = a}, that yields the ideal prior on the d-step SD.
#' @export
prior_conversion <- function(d, prior, p) {
  Cp <- (d^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))
  prior_q <- list(a = prior$a, u = (prior$u * (1 / sqrt(Cp))))
  prior_q
}


dgTMatrix_wrapper <- function(matrix) {
  result <- as(as(as(matrix, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  result
}



#' Roxygen commands
#'
#' @useDynLib OSplines
#'
dummy <- function() {
  return(NULL)
}
