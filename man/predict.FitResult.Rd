% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_post_fit.R
\name{predict.FitResult}
\alias{predict.FitResult}
\title{To predict the GP component in the fitted model, at the locations specified in `newdata`.}
\usage{
\method{predict}{FitResult}(
  object,
  newdata = NULL,
  variable,
  degree = 0,
  include.intercept = TRUE,
  only.samples = FALSE
)
}
\arguments{
\item{object}{The fitted object from the function `model_fit`.}

\item{newdata}{The dataset that contains the locations to be predicted for the specified GP. Its column names must include `variable`.}

\item{variable}{The name of the variable to be predicted, should be in the `newdata`.}

\item{degree}{The degree of derivative that the user specifies for inference. Only applicable for a GP in the `IWP` type.}

\item{include.intercept}{A logical variable specifying whether the intercept should be accounted when doing the prediction. The default is TRUE. For Coxph model, this 
variable will be forced to FALSE.}

\item{only.samples}{A logical variable indicating whether only the posterior samples are required. The default is FALSE, and the summary of posterior samples will be reported.}
}
\description{
To predict the GP component in the fitted model, at the locations specified in `newdata`.
}
