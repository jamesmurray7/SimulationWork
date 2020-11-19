#' ###
#' Simulating joint data to fit a joint model to.
#' This only uses a random intercept for the L.A.
#' and fairly trivial covariate(s), until it works!
#' ###

# Prerequisites ------------------------------------------------------------
dev.off()
rm(list = ls())
library(MASS)
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(survival)

# Setting-out the scenario ------------------------------------------------
# Some trial, higher value of outcome is worse.
# Binary covariate is receiving treatment
# Six treatment times (t)

# Single run --------------------------------------------------------------
# (will functionise afterwards)
b0 <- 40 # Overall mean (intercept)
b1 <- -10 # Perform 10 units better
sigma.i <- 1.5 # Subject R.E (intercept \pm this) (SD)
sigma.e <- 2.5 # Measurement error (SD)
t <- 0:5 # Treatment times
m <- 500 # Number of subjects
n_i <- 6 # Number of measurements (=length(t))
N <- m*n_i

# Data
x <- rep(rbinom(m, 1, 0.5), each = n_i)
U <- rep(rnorm(m, 0, sigma.i), each = n_i) # Random effects
tt <- rep(t, m)
epsilon <- rnorm(N, 0, sigma.e)

Y <- b0 + U + b1 * x + epsilon
long_dat <- data.frame(x, tt, Y)
