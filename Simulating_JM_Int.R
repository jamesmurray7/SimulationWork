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
length(which(survtime > 5))/length(survtime) # ~ 50% experience event with lambda = 1

ratec <- 0.05
censor <- rexp(m, ratec)
id <- 1:m
surv_dat <- data.frame(id, trt, survtime, censor)

# Termination and censoring time
surv_dat$time <- pmin(survtime, censor, max(t)) # When does the profile stop?
surv_dat$status <- ifelse(surv_dat$censor < surv_dat$survtime, 0, 1) # Status (1:died)
surv_dat <- surv_dat[, c("id", "trt", "time", "status")]

summary(coxph(Surv(time, status) ~ trt, data = surv_dat))

# Longitudinal part and survival part are therefore
long_dat %>% head(10)
surv_dat %>% head(10)

joint_dat <- left_join(long_dat, surv_dat, by = "id")

# Cast to class "jointdata"


joineR_joint_dat <- joineR::jointdata(
  longitudinal = long_dat,
  survival = surv_dat,
  baseline = surv_dat[, c("id", "trt")],
  id.col = "id",
  time.col = "tt"
) 

fit <- joineR::joint(joineR_joint_dat, 
              long.formula = Y ~ trt,
              surv.formula = Surv(time, status) ~ trt,
              model = "int",
              sepassoc = F, max.it = 50,
              verbose = T)
summary(fit) # 

