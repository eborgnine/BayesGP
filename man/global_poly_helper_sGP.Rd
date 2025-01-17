% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_utility.R
\name{global_poly_helper_sGP}
\alias{global_poly_helper_sGP}
\title{Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix) of sGP}
\usage{
global_poly_helper_sGP(refined_x, a, m, initial_location = NULL)
}
\arguments{
\item{refined_x}{A vector of locations to evaluate the sB basis}

\item{a}{The frequency of sGP.}

\item{m}{The number of harmonics to consider}
}
\value{
A matrix with i,j componet being the value of jth basis function
value at ith element of x, the ncol should equal to (2*m), and nrow
should equal to the number of elements in x
}
\description{
Constructing and evaluating the global polynomials, to account for boundary conditions (design matrix) of sGP
}
