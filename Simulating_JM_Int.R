#' ###
#' Simulating joint data to fit a joint model to.
#' This only uses a random intercept for the L.A.
#' and fairly trivial covariate(s), until it works!
#' Taking exponential baseline hazard
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
trt <- rbinom(m, 1, 0.5) 
U_int <- rnorm(m, 0, sigma.i)

id <- rep(1:m, each = n_i)
x1 <- rep(trt, each = n_i) # 'extending down' for all times
U <- rep(U_int, each = n_i) # Random effects
tt <- rep(t, m)
epsilon <- rnorm(N, 0, sigma.e)

# Longitudinal part ----

Y <- b0 + U + b1 * x1 + epsilon
long_dat <- data.frame(id, x1, tt, Y)

summary(lmer(Y ~ x1 + (1|id), data = long_dat))

# Survival part ----
lambda <- 1 # BL Hazard

uu <- runif(m)
survtime <- -log(uu)/(lambda * exp(b1*trt))
length(which(survtime > 6))/length(survtime)
