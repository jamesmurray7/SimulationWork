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
lambda <- 0.005
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

long_data <- data.frame(id = rep(id, each = n_i), time, Xl, Y)

summary(lmer(Y ~ x1l + x2l + x3l + time + (1|id), data = long_data)) # Cool!

# Survival part //
Xs <- model.matrix(~x1+x3-1) # Only considering binary and continuous
Bs <- matrix(c(b1_s, b3_s), nrow = 1)

uu <- runif(m)
tt <- -log(uu)/(lambda * exp(Xs %*% t(Bs) + U_int)) 
length(which(tt > max(t)))/length(tt) # % who experience event

# Censoring and truncation
censor <- rexp(m, 0.001)
tau <- max(time)
survtime <- pmin(tt, censor, tau) # time to output
status <- ifelse(survtime == tt, 1, 0)

surv_data <- data.frame(id, x1, x3, survtime, status)

summary(coxph(Surv(survtime, status) ~ x1 + x3, data = surv_data))

# Single-run joint model ----
long_data %>% head(10)
surv_data %>% head(10)

temp <- left_join(long_data, surv_data, "id"); head(temp)

long_data2 <- temp %>% 
  filter(time <= survtime) %>% 
  dplyr::select(names(long_data))

# Cast to class "jointdata"

jd <- joineR::jointdata(
  longitudinal = long_data2,
  survival = surv_data,
  baseline = surv_data[, c("id", "x1", "x3")],
  id.col = "id",
  time.col = "time"
) 

fit <- joineR::joint(jd, 
              long.formula = Y ~ x1l + x2l +x3l + time,
              surv.formula = Surv(survtime, status) ~ x1 + x3,
              model = "int")
summary(fit) 


# Functionise -------------------------------------------------------------
# Just random intercept again!

joint_sim <- function(m = 250, n_i = 5, 
                      Bl = c(40, -10, 5, 15, 0.1), # Longit: Intercept, binary, factor2-3, continuous
                      Bs = c(-0.3, 0.05), # Survival: log-odds binary and continuous,
                      sigma.i = 1.5, sigma.e = 2.5,
                      lambda = 0.005){
  # Common parameters //
  N <- m * n_i
  id <- 1:m
  tau <- n_i-1
  # SDs //
  sigma.e <- 2.5 # measurement
  sigma.i <- 1.5 # intercept
  # Error terms //
  U_int <- rnorm(m, 0, sigma.i)
  epsilon <- rnorm(N, 0, sigma.e)
  # Longitudinal coefficients //
  Bl <- matrix(Bl, nrow = 1) # Cast to matrix
  # Survival coefficients //
  lambda <- 0.005
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
  
  long_data <- data.frame(id = rep(id, each = n_i), time, Xl, Y)
  
  summary(lmer(Y ~ x1l + x2l + x3l + time + (1|id), data = long_data)) # Cool!
  
  # Survival part //
  Xs <- model.matrix(~x1+x3-1) # Only considering binary and continuous
  Bs <- matrix(c(b1_s, b3_s), nrow = 1)
  
  uu <- runif(m)
  tt <- -log(uu)/(lambda * exp(Xs %*% t(Bs) + U_int)) 
  length(which(tt > max(t)))/length(tt) # % who experience event
  
  # Censoring and truncation
  censor <- rexp(m, 0.001)
  tau <- max(time)
  survtime <- pmin(tt, censor, tau) # time to output
  status <- ifelse(survtime == tt, 1, 0)
  
  surv_data <- data.frame(id, x1, x3, survtime, status)
  
  # Extra output - number of events
  pc_events <- length(which(survtime < tau))/m * 100
  
  return(list(long_dat, surv_dat, pc_events))
  
}

temp <- joint_sim()
summary(lmer(Y ~ x1l + x2l + x3l + time + (1|id), data = temp[[1]]))
summary(coxph(Surv(survtime, status) ~ x1 + x3, data = temp[[2]]))


# Separate investigation --------------------------------------------------
# Should illustrate need for JM
separate_fits <- function(df){
  lmm_fit <- lmer(Y ~ x1l + x2l + x3l + time + (1|id), data = df[[1]])
  surv_fit <- coxph(Surv(survtime, status) ~ x1 + x3, data = df[[2]])
  return(
    list(lmm_fit, surv_fit)
  )
}

pb <- progress::progress_bar$new(total = 1000)
longit_beta <- data.frame(beta0 = NA, beta1 = NA, beta22 = NA, beta23 = NA, beta3 = NA, sigma.e = NA, sigma.u = NA)
surv_beta <- data.frame(beta1s = NA, beta3s = NA)
pc_events <- c()

for(i in 1:1000){
  dat <- joint_sim()
  pc_events[i] <- dat[[3]]
  fits <- separate_fits(dat)
  long_coefs <- fits[[1]]@beta[1:5]
  long_sigma.e <- sigma(fits[[1]])
  long_sigma.u <- as.numeric(attr(VarCorr(fits[[1]])$id, "stddev"))
  longit_beta[i,] <- c(long_coefs, long_sigma.e, long_sigma.u)
  surv_beta[i, ] <- as.numeric(fits[[2]]$coefficients)
  pb$tick()
}


ex <- expression
to_plot <- cbind(longit_beta, surv_beta, pc_events) %>% tibble %>% 
  gather("parameter", "estimate") %>% 
  mutate(param = factor(parameter, levels = c("beta0", "beta1", "beta22", "beta23", "beta3", "sigma.e", "sigma.u",
                                              "beta1s", "beta3s", "pc_events"),
                        labels = c(ex(beta[0]), ex(beta[1]), ex(beta[22]), ex(beta[23]), ex(beta[3]),
                                   ex(sigma[e]), ex(sigma[u]), ex(beta[1*"S"]), ex(beta[3*"S"]), ex("Events")))
         )

plot_lines <- to_plot %>% distinct(param)
plot_lines$xint <- c(40, -10, 5, 15, 0.1, 2.5, 1.5, -0.3, 0.05, NA)

to_plot %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .2) + 
  geom_vline(data = plot_lines, aes(xintercept = xint), colour = "blue", alpha = .5, lty = 3) + 
  facet_wrap(~param, scales = "free", nrow = 5, ncol = 2, labeller = label_parsed) + 
  labs(title = "Separate investigation", x = "Estimate")
ggsave("./JM-sims-plots/Separate_Investigation.png")



# Joint investigation -----------------------------------------------------
library(joineR)

long_dat <- joint_sim()[[1]]
surv_dat <- joint_sim()[[2]]

# Single-run
temp <- left_join(long_dat, surv_dat, "id")

long_dat2 <- temp %>% 
  filter(time <= survtime) %>% 
  dplyr::select(names(long_dat))

jd <- jointdata(
  longitudinal = long_dat2,
  survival = surv_dat,
  id.col = "id",
  time.col = "time",
  baseline = surv_dat[,c("id", "x1", "x3")]
)

joint_fit <- joint(jd,
      long.formula = Y ~ x1l + x2l + x3l + time,
      surv.formula = Surv(survtime, status) ~ x1 + x3,
      model = "int") # Sepassoc doesn't matter as only one L.A.

summary(joint_fit)

# Function - allowing us to change sample sizes
joint_fit <- function(m, n_i){
  dat <- joint_sim(m = m, n = n_i)
  long_data <- left_join(dat[[1]], dat[[2]], "id") %>% 
    filter(time <= survtime) %>% 
    select(names(dat[[1]]))
  jd <- jointdata(
    longitudinal = long_data,
    survival = dat[[2]],
    time.col = "time", id.col = "id",
    baseline = dat[[2]][, c("id", "x1", "x3")]
  )
  fit <- joint(jd,
               long.formula = Y ~ x1l + x2l + x3l + time,
               surv.formula = Surv(survtime, status) ~ x1 + x3,
               model = "int")
  return(fit)
}

smallfn <- function(x){
  replicate(x, joint_fit(100, 5), simplify = F)
}
medfn <- function(x){
  replicate(x, joint_fit(250, 6), simplify = F)
}
largefn <- function(x){
  replicate(x, joint_fit(500, 5), simplify = F)
}

library(furrr)
plan(multicore)

small_results <- future_map(500, smallfn)
plan(sequential)
plan(multicore)
med_results <- future_map(500, medfn)
plan(sequential)
plan(multicore)
large_results <- future_map(500, largefn)
plan(sequential)

# Extract joint fits and compare.
extract_coefs <- function(fit){
  converged <- fit$convergence
  # Longitudinal parameters
  longits <- t(as.numeric(t(fit$coefficients$fixed$longitudinal))[1:5])
  sigma.e <- sqrt(fit$sigma.z)
  sigma.u <- as.numeric(sqrt(fit$sigma.u))
  # Survival
  surv <- t(as.numeric(fit$coefficients$fixed$survival))
  # Latent association, gamma
  latent_association <- as.numeric(fit$coefficients$latent)
  
  return(data.frame(converged, longits, sigma.e, sigma.u, surv, latent_association))
}

small_results2 <- tibble(small_results[[1]]) %>% mutate(params = map(`small_results[[1]]`, extract_coefs))
med_results2 <- tibble(med_results[[1]]) %>% mutate(params = map(`med_results[[1]]`, extract_coefs))
large_results2 <- tibble(large_results[[1]]) %>% mutate(params = map(`large_results[[1]]`, extract_coefs))

ex <- expression
param_plot <- function(df){
  to_plot <- df %>% 
    unnest(params) %>% 
    rename(b0 = X1, b1 = X2, b22 = X3, b23 = X4, b3 = X5, b1s = `X1.1`, b3s = `X2.1`) %>% 
    select(-1) %>% 
    gather("parameter", "Estimate", -converged) %>% 
    filter(converged) %>% 
    mutate(
      param = factor(parameter,
                     levels = c("b0", "b1", "b22", "b23", "b3", "sigma.u", "sigma.e",
                                "b1s", "b3s", "latent_association"),
                     labels = c(ex(beta[0]), ex(beta[1]), ex(beta[22]),
                                ex(beta[23]), ex(beta[3]), ex(sigma[u]), ex(sigma[e]),
                                ex(beta[1*"S"]), ex(beta[3*"S"]), ex(gamma[0])))
    )
  
  xints <- to_plot %>% distinct(param)
  xints$xint <- c(40, -10, 5, 15, 0.1, 2.5, 1.5, -0.3, 0.05, 1)

  to_plot %>% 
    ggplot(aes(x = Estimate)) + 
    geom_density(alpha = .2, fill = "grey20") + 
    geom_vline(data = xints, aes(xintercept = xint), colour = "blue", lty = 1, alpha = .8) + 
    facet_wrap(~param, scales = "free", ncol = 2, nrow = 5, labeller = label_parsed) + 
    theme(strip.background = element_blank(),
          strip.text.x = element_text(colour = "black", size = 12))
}

param_plot(small_results2)
ggsave("./JM-sims-plots/SmallSample.png")
param_plot(med_results2)
ggsave("./JM-sims-plots/MediumSample.png")
param_plot(large_results2)
ggsave("./JM-sims-plots/LargeSample.png")
