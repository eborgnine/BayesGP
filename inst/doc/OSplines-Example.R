## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(OSplines)
library(tidyverse)
library(Matrix)
library(TMB)
library(aghq)

## -----------------------------------------------------------------------------
data <- INLA::Munich %>% select(rent, year)
head(data, n = 5)

## -----------------------------------------------------------------------------
### Initialization
data$t <- data$year - min(data$year)
length(unique(data$t))

## -----------------------------------------------------------------------------
### Use the observed locations as knots
knots_observed <- sort(unique(data$t))

### Use equally spaced knots, with 4 years as spacing
knots_equal <- sort(seq(from = min(data$t), to = max(data$t), by = 4))


## -----------------------------------------------------------------------------
knots_observed[1]
knots_equal[1]

## -----------------------------------------------------------------------------
B1 <- local_poly(knots = knots_observed, 
                 refined_x = data$t,
                 p = 2)
B2 <- local_poly(knots = knots_equal, 
                 refined_x = data$t,
                 p = 2)

## -----------------------------------------------------------------------------
P1 <- compute_weights_precision(x = knots_observed)
P2 <- compute_weights_precision(x = knots_equal)

## -----------------------------------------------------------------------------
dim(B1)
dim(P1)
dim(B2)
dim(P2)

## -----------------------------------------------------------------------------
X <- global_poly(x = data$t, p = 2)
dim(X)

## -----------------------------------------------------------------------------
tmbdat1 <- list(
          # Design matrix
          X = as(X,'dgTMatrix'),
          B = as(B1,'dgTMatrix'),
          P = as(P1,'dgTMatrix'),
          logPdet = as.numeric(determinant(P1,logarithm = T)$modulus),
          # Response
          y = data$rent,
          # PC Prior params
          # (u1, alpha1) for sigma_s
          # (u2, alpha2) for sigma
          u1 = 1,
          alpha1 = 0.5,
          u2 = 1,
          alpha2 = 0.5,
          betaprec = 0.01
        )

tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B1)))), # W = c(U,beta); U = Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )

ff <- TMB::MakeADFun(
          data = tmbdat1,
          parameters = tmbparams,
          random = "W",
          DLL = "OSplines",
          silent = TRUE
        )

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
mod1 <- aghq::marginal_laplace_tmb(ff,4,c(0,0))


## -----------------------------------------------------------------------------
tmbdat2 <- list(
          # Design matrix
          X = as(X,'dgTMatrix'),
          B = as(B2,'dgTMatrix'),
          P = as(P2,'dgTMatrix'),
          logPdet = as.numeric(determinant(P1,logarithm = T)$modulus),
          # Response
          y = data$rent,
          # PC Prior params
          # (u1, alpha1) for sigma_s
          # (u2, alpha2) for sigma
          u1 = 1,
          alpha1 = 0.5,
          u2 = 1,
          alpha2 = 0.5,
          betaprec = 0.01
        )

tmbparams <- list(
          W = c(rep(0, (ncol(X) + ncol(B2)))), # W = c(U,beta); U = Spline coefficients
          theta1 = 0, # -2log(sigma)
          theta2 = 0
        )

ff <- TMB::MakeADFun(
          data = tmbdat2,
          parameters = tmbparams,
          random = "W",
          DLL = "OSplines",
          silent = TRUE
        )

# Hessian not implemented for RE models
ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
mod2 <- aghq::marginal_laplace_tmb(ff,4,c(0,0))

## -----------------------------------------------------------------------------
samps1 <- sample_marginal(mod1, M = 3000)
global_samps1 <- samps1$samps[(ncol(B1) + 1):nrow(samps1$samps),]
coefsamps1 <- samps1$samps[1:ncol(B1),]

samps2 <- sample_marginal(mod2, M = 3000)
global_samps2 <- samps2$samps[(ncol(B2) + 1):nrow(samps2$samps),]
coefsamps2 <- samps2$samps[1:ncol(B2),]

## -----------------------------------------------------------------------------
f1 <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                       knots = knots_observed, 
                       refined_x = seq(from = min(data$t), to = max(data$t), by = 1),
                       p = 2, degree = 0)

f2 <- compute_post_fun(samps = coefsamps2, global_samps = global_samps2, 
                       knots = knots_equal, 
                       refined_x = seq(from = min(data$t), to = max(data$t), by = 1),
                       p = 2, degree = 0)


f1pos <- extract_mean_interval_given_samps(f1)
f2pos <- extract_mean_interval_given_samps(f2)

## -----------------------------------------------------------------------------
f1pos %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), color = "blue") + 
  geom_ribbon(aes(ymin = plower, ymax = pupper), fill = "orange", alpha = 0.3) +
  geom_point(aes(x = t, y = rent), alpha = 0.1, data = data) +
  theme_classic() + xlab("year since 1981") + ylab("rent") +
  ggtitle("Using observed knots") + ylim(c(5,12))

f2pos %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), color = "blue") + 
  geom_ribbon(aes(ymin = plower, ymax = pupper), fill = "orange", alpha = 0.3) +
  geom_point(aes(x = t, y = rent), alpha = 0.1, data = data) +
  theme_classic() + xlab("year since 1981") + ylab("rent") +
  ggtitle("Using equally spaced knots") + ylim(c(5,12))

## -----------------------------------------------------------------------------
f1deriv <- compute_post_fun(samps = coefsamps1, global_samps = global_samps1, 
                       knots = knots_observed, 
                       refined_x = seq(from = min(data$t), to = max(data$t), by = 1),
                       p = 2, degree = 1)

f2deriv <- compute_post_fun(samps = coefsamps2, global_samps = global_samps2, 
                       knots = knots_equal, 
                       refined_x = seq(from = min(data$t), to = max(data$t), by = 1),
                       p = 2, degree = 1)


f1derivpos <- extract_mean_interval_given_samps(f1deriv)
f2derivpos <- extract_mean_interval_given_samps(f2deriv)

## -----------------------------------------------------------------------------
f1derivpos %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), color = "blue") + 
  geom_ribbon(aes(ymin = plower, ymax = pupper), fill = "orange", alpha = 0.3) +
  theme_classic() + xlab("year since 1981") + ylab("rent: rate of change") +
  ggtitle("Using observed knots")

f2derivpos %>% ggplot(aes(x = x)) + geom_line(aes(y = mean), color = "blue") + 
  geom_ribbon(aes(ymin = plower, ymax = pupper), fill = "orange", alpha = 0.3) +
  theme_classic() + xlab("year since 1981") + ylab("rent: rate of change") +
  ggtitle("Using equally spaced knots")

