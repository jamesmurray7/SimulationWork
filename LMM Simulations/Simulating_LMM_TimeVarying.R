#' ### 
#' Updating LMM sim to include a time varying covariate in two different ways
#' (Unsure which is correct, though!)
#' ###


# Prerequisites -----------------------------------------------------------
rm(list = ls())
dev.off()
library(lme4)
library(MASS)
library(dplyr)

#' ###
#' Scenario:
#' Some trial (higher value = worse)
#' Treatment (x1); age (x2) at baseline; blood pressure at different times (x3(t))
#' ###

# Random intercept: Method 1 ----------------------------------------------
# Global parameters //
m <- 500
n_i <- 6
N <- m * n_i
time <- 0:(n_i-1)
# SDs //
sigma.e <- 2
sigma.u <- 3.5
# Coefficients //
b0 <- 40
b1 <- -5
b2 <- 0.5
b3 <- -0.5
# Covariates //
x1 <- rep(rbinom(m, 1, 0.5), each = n_i)
x2 <- rep(floor(rnorm(m, 70, 10)), each = n_i)
x3 <- rnorm(N, 60, 20)
# Error terms //
U <- rep(rnorm(m, 0, sigma.u), each = n_i)
epsilon <- rnorm(N, 0, sigma.e)
# Outcome, data and model //
Y <- b0 + U + b1 * x1 + b2 * x2 + b3 * x3 + epsilon

long_data <- data.frame(
  id = rep(1:m, each = n_i),
  time = rep(time, m),
  x1, x2, x3, Y
)

summary(lmer(Y ~ x1 + x2 + x3 + time + (1|id), data = long_data))

# Random intercept: Method 2 ----------------------------------------------
rm(x3)
# All other covariates staying the same, now implementing structure to x3 //
grad <- -2
x30 <- rnorm(m, 60, 15)
# Initialise vector //
x3 <- vector("numeric", length = N)
x3 <- rep(x30, each = n_i)
# Time component //
tt <- rep(time, m)
# Rewrite X3 as a function of time //
x3 <- x3 + tt * rnorm(N, grad, 1)
# Outcome, data and model //
Y <- b0 + U + b1 * x1 + b2 * x2 + b3 * x3 + epsilon

long_data <- data.frame(
  id = rep(1:m, each = n_i),
  time = rep(time, m),
  x1, x2, x3, Y
)

summary(lmer(Y ~ x1 + x2 + x3 + time + (1|id), data = long_data))


# Random Slope  -----------------------------------------------------------

# Will be much the same as non-time-varying - do at later date (copy paste exercise!)
