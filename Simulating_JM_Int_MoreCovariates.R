#' ###
#' Simulating joint data to fit a joint model to.
#' This only uses a random intercept for the L.A.
#' Expanding upon Simulating_JM_Int.R by introducing more covariates
#' Taking exponential, and later Weibull baseline hazard
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
# Some diabetes trial, higher value of outcome is worse.
# Binary covariate is receiving treatment (yes = good)
# Factor covariate is BMI category at baseline (fat = bad)
# Continuous covariate is age at baseline (older = bad)
# Six treatment times (t)

# Single run --------------------------------------------------------------
# (will functionise afterwards)
# Common parameters //
t <- 6
m <- 500
n_i <- 6
N <- m * n_i
# SDs //
sigma.e <- 2.5 # measurement
sigma.i <- 1.5 # intercept
# Error terms //
U_int <- rnorm(m, 0, sigma.i)
epsilon <- rnorm(N, 0, sigma.e)
# Longitudinal coefficients //
b0 <- 40
b1_l <- -10
b22_l <- 5
b23_l <- 15
b3_l <- 0.1
Bl <- matrix(c(b0, b1_l, b22_l, b23_l, b3_l), nrow = 1) # Cast to matrix
# Survival coefficients //
lambda <- 0.05
b1_s <- -0.3 # log-odds associated with having treatment (30% HR reduction)
b3_s <- 0.05 # log-odds associated with one unit increase age (5% HR increase)
# Data - baseline // 
x1 <- rbinom(m, 1, 0.5) # Treatment received
x2 <- gl(3, 1, m) # Factor
x3 <- floor(rnorm(m, 65, 7)) # Age
id <- 1:m

# Longitudinal part //
x1l <- rep(x1, each = n_i)
x2l <- rep(x2, each = n_i)
x3l <- rep(x3, each = n_i)
Xl <- model.matrix(~x1l+x2l+x3l)
Ul <- rep(U_int, each = n_i)
time <- rep(0:(t-1), m)

Y <- Xl %*% t(Bl) + Ul + epsilon

long_data <- data.frame(id = rep(id, each = n_i), time, x1l, x2l, x3l, Y)

summary(lmer(Y ~ x1l + x2l + x3l + time + (1|id), data = long_data)) # Cool!

# Survival part ----
lambda <- 0.01 # BL Hazard

uu <- runif(m)
survtime <- -log(uu)/(lambda * exp(b1*trt + U_int)) # This is wrong - needs its own coeff - fixed in function
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
summary(fit) 


# Functionise -------------------------------------------------------------
# Just random intercept and one binary covariate, again.

joint_sim <- function(m = 500, n_i = 6, 
                      b0 = 40, b1 = -10, b1s = -0.5,
                      sigma.i = 1.5, sigma.e = 2.5,
                      lambda = 0.05){
  # Set out variables
  N <-  m * n_i
  id <- 1:m
  time <- 0:(n_i-1)
  tau <- max(time) 
  U_int <- rnorm(m, 0, sigma.i) # Random effects
  x <- rbinom(m, 1, 0.5) # Treatment assignment per id
  
  # Longitudinal part //
  
  xl <- rep(x, each = n_i)
  U <- rep(U_int, each = n_i)
  epsilon <- rnorm(N, 0, sigma.e)
  Y <- b0 + U + b1 * xl + epsilon
  
  long_dat <- data.frame(id = rep(id, each = n_i),
                         time = rep(time, m),
                         xl, Y)
  
  # Survival part     //
  
  # Survival times
  u <- runif(m)
  
  tt <- -log(u)/(lambda * exp(b1s * x))
  
  # Censoring and truncation
  rateC <- 0.05
  censor <- rexp(m, rateC)
  survtime <- pmin(tt, censor, tau) # time to output
  status <- ifelse(survtime == tt, 1, 0)
  
  surv_dat <- data.frame(id, x, survtime, status)
  
  # Extra output - number of events
  pc_events <- length(which(survtime < tau))/m * 100
  
  return(list(long_dat, surv_dat, pc_events))
  
}

temp <- joint_sim()
lmer(Y ~ xl + time + (1|id), data = temp[[1]])
coxph(Surv(survtime, status) ~ x, data = temp[[2]]) # Look


# Investigate -------------------------------------------------------------

# Separate investigation ----

separate_fits <- function(df){
  lmm_fit <- lmer(Y ~ xl + time + (1|id), data = df[[1]])
  surv_fit <- coxph(Surv(survtime, status) ~ x, data = df[[2]])
  return(
     list(lmm_fit, surv_fit)
  )
}

pb <- progress::progress_bar$new(total = 1000)
longit_beta <- surv_beta <- pc_events <-  c()
for(i in 1:1000){
  dat <- joint_sim()
  pc_events[i] <- dat[[3]]
  fits <- separate_fits(dat)
  longit_beta[i] <- fits[[1]]@beta[2]
  surv_beta[i] <- fits[[2]]$coefficients
  pb$tick()
}

mean(longit_beta); mean(surv_beta)

# Confirm it's got the right things in this separate case.
data.frame(lmm = longit_beta, cox = surv_beta, pc_events) %>%
  gather("outcome", "value") %>% 
  ggplot(aes(x = value)) +
  geom_density(fill = "grey20", alpha = .5) + 
  facet_wrap(~outcome, scales = "free")

# Joint investigation ----
library(joineR)

long_dat <- joint_sim()[[1]]
surv_dat <- joint_sim()[[2]]

# Remove where IDs have failed

temp <- left_join(long_dat, surv_dat, "id")

long_dat2 <- temp %>% 
  filter(time <= survtime) %>% 
  dplyr::select(names(long_dat))

jd <- jointdata(
  longitudinal = long_dat2,
  survival = surv_dat,
  id.col = "id",
  time.col = "time",
  baseline = surv_dat[,c("id", "x")]
)

joint_fit <- joint(jd,
      long.formula = Y ~ xl + time,
      surv.formula = Surv(survtime, status) ~ x,
      model = "int") # Sepassoc doesn't matter as only one L.A.

summary(joint_fit)


#' ###
#' TO DO:
#' Functionise the above joint fit and 
#' parse through the list's outputs to obtain estimates
#' and do some plots. After, move on to more complex 
#' scenario (3x covariates) and then random slope
#' ###

# Functionise joint fit ---------------------------------------------------

joint_fit <- function(...){
  # Initialise long and survival parts
  dat <- joint_sim(...)
  temp <- left_join(dat[[1]], dat[[2]], "id")
  long <- temp %>% filter(time <= survtime) %>% dplyr::select(names(dat[[1]]))
  surv <- dat[[2]]
  # Cast to class joint data //
  jd <- jointdata(
    longitudinal = long, survival = surv,
    id.col = "id", time.col = "time",
    baseline = surv[, c("id", "x")]
  )
  # Fit joint model //
  fit <- joint(jd,
               long.formula = Y ~ xl + time,
               surv.formula = Surv(survtime, status) ~ x,
               model = "int")
  # Extract parameters of interest //
  epsilon <- sqrt(fit$sigma.z) # Random error SD
  U <- sqrt(as.numeric(fit$sigma.u)) # Random effects SD
  beta_l <- fit$coefficients$fixed$longitudinal[2,1]
  beta_s <- fit$coefficients$fixed$survival
  # Return data frame of these coefficients
  return(
    data.frame(beta_l, beta_s, U, epsilon)
  )
}

fits <- replicate(100, joint_fit(), simplify = F) # Default 500 x 6 data
fits_smallsample <- replicate(100, joint_fit(m = 100, n_i = 5), simplify = F)
fits_largersample <- replicate(100, joint_fit(m = 750, n_i = 10), simplify = F)

# Plot these lists of model fits
plots <- list()
plotfn <- function(fitlist){
  plot.out <- fitlist %>% 
    bind_rows %>% 
    gather("parameter", "estimate") %>% 
    mutate(
      param = factor(parameter, levels = c("beta_l", "beta_s", "epsilon", "U"),
                     labels = c(expression(beta[longit]), expression(beta[surv]),
                                expression(sigma[e]), expression(sigma[U])))
    ) %>% 
    ggplot(aes(x = estimate)) + 
    geom_density(colour = "grey20", alpha = .2) + 
    facet_wrap(~param, ncol = 4, nrow = 1, scales = "free", labeller = label_parsed) + 
    theme(
      strip.text = element_text(size = 12, colour = "black"),
      strip.background = element_blank()
    )
  return(plot.out)
}

ggpubr::ggarrange(
  plotfn(fits_smallsample),
  plotfn(fits),
  plotfn(fits_largersample),
  ncol = 1
)
ggsave("./JM-sims-plots/SampleSizeInt.png")
  