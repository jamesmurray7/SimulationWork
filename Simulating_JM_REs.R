#' ###
#' Simulating joint data to fit a joint model to.
#' Taking work started in Simulating_JM_Int_MoreCovariates.R
#' and adding random slope to model specifications.
#' Latent association therefore dependent on R.E: Intercept and Slope (gamma_1 and gamma_2)
#' in Henderson 2000.
#' ###

# Prerequisites ------------------------------------------------------------
dev.off()
rm(list = ls())
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

# 'Global' parameters //
n_i <- 6 
m <- 250
N <- n_i * m
rho <- 0.3
sigma.i <- 3   # Random intercept SD
sigma.s <- 2   # Random slope SD
sigma.e <- 1.5 # Epsilon SD

# Covariance matrix //
sigma <- matrix(
  c(sigma.i ^ 2, rho * sigma.i * sigma.s,
  rho * sigma.i * sigma.s, sigma.s^2), nrow = 2, byrow = T
)

# Generate Random effects //
RE <- MASS::mvrnorm(m, c(0,0), sigma)
Ui <- RE[, 1]; Us <- RE[, 2] # Intercept; Slope

# Covariates - baseline //
id <- 1:m
x1 <- rbinom(m, 1, 0.5) # Treatment received
x2 <- gl(3, 1, m) # Factor
x3 <- floor(rnorm(m, 65, 10)) # Age

# Longitudinal part //
# Coefficients
b0 <- 40; b1 <- -10; b22 <- 5; b23 <- 15; b3 <- 0.1
Bl <- matrix(c(b0, b1, b22, b23, b3), nrow = 1) # Cast to matrix
# Baseline covariates
x1_l <- rep(x1, each = n_i)
x2_l <- rep(x2, each = n_i)
x3_l <- rep(x3, each = n_i)
Xl <- model.matrix(~x1_l+x2_l+x3_l)
time <- rep(0:(n_i-1), m)
# REs
U1l <- rep(Ui, each = n_i)
U2l <- rep(Us, each = n_i)
epsilon <- rnorm(N, 0, sigma.e)
# Response
Y <- Xl %*% t(Bl) + U1l + U2l * time + epsilon
# Data and quick model
long_data <- data.frame(id = rep(id, each = n_i), x1_l, x2_l, x3_l, time, Y)
summary(lmer(Y ~ x1_l + x2_l + x3_l + time + (1+time|id), data = long_data)) # Cool!

# Survival part //
lambda <- 0.05
b1s <- -0.3 # log-odds associated with having treatment (30% HR reduction)
b3s <- 0.05 # log-odds associated with one unit increase age (5% HR increase)  
Bs <- matrix(c(b1s, b3s), nrow = 1)
Xs <- model.matrix(~ x1 + x3 - 1)
  
# Simulate survival times
uu <- runif(m)
tt <- -log(uu)/(lambda * exp(Xs %*% t(Bs) + Ui + Us))
# Censoring
censor <- rexp(m, 0.01)
survtime <- pmin(tt, censor, 5)
status <- ifelse(survtime == tt, 1, 0)

surv_data <- data.frame(id, x1, x3, survtime, status)

summary(coxph(Surv(survtime, status) ~ x1 + x3, data = surv_data)) # Way further off than just R.I!

# Cast to class "jointdata"

jd <- joineR::jointdata(
  longitudinal = long_data,
  survival = surv_data,
  time.col = "time",
  id.col = "id",
  baseline = surv_data[, c("id", "x1", "x3")]
)

summary(joineR::joint(
  data = jd,
  long.formula = Y ~ x1_l + x2_l + x3_l + time,
  surv.formula = Surv(survtime, status) ~ x1 + x3
)) # Cool cool!


# Functionise -------------------------------------------------------------
# Just random intercept again!

joint_sim <- function(m = 200, n_i = 5, 
                      Bl = c(40, -10, 5, 15, 0.1), # Longit: Intercept, binary, factor2-3, continuous
                      Bs = c(-0.3, 0.05), # Survival: log-odds binary and continuous,
                      sigma.i = 3, sigma.s = 2, rho = 0.3, # RE parameters
                      sigma.e = 1.5, # Error parameter
                      lambda = 0.05){
  # 'Global' parameters //
  N <- n_i * m
  id <- 1:m
  tau <- n_i-1 # UL time variable.
  rho <- 0.3
  sigma.i <- 3   # Random intercept SD
  sigma.s <- 2   # Random slope SD
  sigma.e <- 1.5 # Epsilon SD
  
  # Covariance matrix //
  sigma <- matrix(
    c(sigma.i ^ 2, rho * sigma.i * sigma.s,
      rho * sigma.i * sigma.s, sigma.s^2), nrow = 2, byrow = T
  )
  
  # Generate Random effects //
  RE <- MASS::mvrnorm(m, c(0,0), sigma)
  Ui <- RE[, 1]; Us <- RE[, 2] # Intercept; Slope
  
  # Covariates - baseline //
  x1 <- rbinom(m, 1, 0.5) # Treatment received
  x2 <- gl(3, 1, m) # Factor
  x3 <- floor(rnorm(m, 65, 10)) # Age
  
  # Longitudinal part //
  # Coefficients
  Bl <- matrix(Bl, nrow = 1) # Cast to matrix
  # Baseline covariates
  x1l <- rep(x1, each = n_i)
  x2l <- rep(x2, each = n_i)
  x3l <- rep(x3, each = n_i)
  Xl <- model.matrix(~x1l+x2l+x3l)
  time <- rep(0:tau, m)
  # REs
  U1l <- rep(Ui, each = n_i)
  U2l <- rep(Us, each = n_i)
  epsilon <- rnorm(N, 0, sigma.e)
  # Response
  Y <- Xl %*% t(Bl) + U1l + U2l * time + epsilon
  # Data and quick model
  long_data <- data.frame(id = rep(id, each = n_i), x1l, x2l, x3l, time, Y)
  
  # Survival part //
  Bs <- matrix(Bs, nrow = 1)
  Xs <- model.matrix(~ x1 + x3 - 1)
  
  # Simulate survival times
  uu <- runif(m)
  tt <- -log(uu)/(lambda * exp(Xs %*% t(Bs) + Ui + Us))
  # Censoring
  censor <- rexp(m, 0.01)
  survtime <- pmin(tt, censor, 5)
  status <- ifelse(survtime == tt, 1, 0)
  
  surv_data <- data.frame(id, x1, x3, survtime, status)
  
  # Extra output - number of events
  pc_events <- length(which(survtime < tau))/m * 100
  
  return(list(long_data, surv_data, pc_events))
  
}

temp <- joint_sim()
summary(lmer(Y ~ x1l + x2l + x3l + time + (1+time|id), data = temp[[1]]))
summary(coxph(Surv(survtime, status) ~ x1 + x3, data = temp[[2]]))


# Separate investigation --------------------------------------------------
# Should illustrate need for JM
separate_fits <- function(df){
  lmm_fit <- lmer(Y ~ x1l + x2l + x3l + time + (1+time|id), data = df[[1]])
  surv_fit <- coxph(Surv(survtime, status) ~ x1 + x3, data = df[[2]])
  return(
    list(lmm_fit, surv_fit)
  )
}

pb <- progress::progress_bar$new(total = 1000)
longit_beta <- data.frame(beta0 = NA, beta1 = NA, beta22 = NA, beta23 = NA, beta3 = NA, 
                          sigma.e = NA, sigma.ui = NA, sigma_us = NA)
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
  mutate(param = factor(parameter, levels = c("beta0", "beta1", "beta22", "beta23", "beta3", 
                                              "sigma.e", "sigma.ui", 'sigma_us',
                                              "beta1s", "beta3s", "pc_events"),
                        labels = c(ex(beta[0]), ex(beta[1]), ex(beta[22]), ex(beta[23]), ex(beta[3]),
                                   ex(sigma[e]), ex(sigma[u*"i"]), ex(sigma[u*"s"]),
                                   ex(beta[1*"S"]), ex(beta[3*"S"]), ex("Percent"*" experiencing"*" events")))
         )

plot_lines <- to_plot %>% distinct(param)
plot_lines$xint <- c(40, -10, 5, 15, 0.1, 1.5, 3, 2, -0.3, 0.05, NA)

to_plot %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .2) + 
  geom_vline(data = plot_lines, aes(xintercept = xint), colour = "blue", alpha = .5, lty = 3) + 
  facet_wrap(~param, scales = "free", nrow = 6, ncol = 2, labeller = label_parsed) + 
  labs(title = "Separate investigation", x = "Estimate")
ggsave("./JM-sims-plots/Separate_Investigation_REs.png")

# Joint investigation -----------------------------------------------------
library(joineR)

long_dat <- joint_sim()[[1]]
surv_dat <- joint_sim()[[2]]

# Single-run ----
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
      sepassoc = T) # Sepassoc doesn't matter as only one L.A.

summary(joint_fit)

# Functionise casting to joint data and fitting joint model ----

simstudy <- tibble(m = c(100, 250, 500), n_i = c(5, 4, 4)) %>% 
  mutate(id = glue::glue("m = {m}, n_i = {n_i}")) %>% 
  group_by(id) %>% 
  nest() %>% 
  ungroup

joint_fit <- function(dat){
  m <- dat$m; n_i <- dat$n_i
  
  pb <- progress::progress_bar$new(total = 250)
  joint_fit_list <- list()
  
  for(i in 1:250){
    df <- joint_sim(m = m, n_i = n_i)
    long_data <- left_join(df[[1]], df[[2]], "id") %>% 
      filter(time <= survtime) %>% 
      select(names(df[[1]]))
    
    jd <- jointdata(
      longitudinal = long_data,
      survival = df[[2]],
      baseline = df[[2]][, c("id", "x1", "x3")],
      id.col = "id", time.col = "time"
    )
    
    joint_fit <- joint(jd,
                       long.formula = Y ~ x1l + x2l + x3l,
                       surv.formula = Surv(survtime, status) ~ x1 + x3,
                       sepassoc = T)
    joint_fit_list[[i]] <- joint_fit
    pb$tick()
  }
  return(joint_fit_list)
}
tic()
simstudy$results <- parallel::mclapply(simstudy$data, joint_fit,
                              mc.preschedule = T, mc.cores = 2) # 53 minutes!
toc()

simstudy$results[[1]][[1]]


extract_parameters <- function(fit){
  Bl <- fit$coefficients$fixed$longitudinal$b1
  sigma.e <- sqrt(fit$sigma.z)
  sigma.u <- sqrt(diag(fit$sigma.u))
  Bs <- as.numeric(fit$coefficients$fixed$survival)
  gammas <- fit$coefficients$latent
  out <- c(Bl, sigma.e, sigma.u, Bs, gammas)
  out <- as.data.frame(t(out))
  names(out) <- c("b0", "b1", "b22", "b23", "b3", "sigma.e", "sigma.ui",
                  "sigma.us", "b1s", "b3s", 'gamma0', "gamma1")
  return(
    out
  )
}
extract_parameters(simstudy$results[[1]][[1]])

small_results <- bind_rows(lapply(simstudy$results[[1]], extract_parameters))
medium_results <- bind_rows(lapply(simstudy$results[[2]], extract_parameters))
large_results <- bind_rows(lapply(simstudy$results[[3]], extract_parameters))

small_results$id <- "m = 100, n_i = 5"
medium_results$id <- "m = 250, n_i = 4"
large_results$id <- "m = 500, n_i = 4"

all_results <- rbind(small_results, medium_results, large_results)

all_results2 <- all_results %>% 
  mutate(id = fct_inorder(id)) %>% 
  gather("Parameter", "Estimate", -id) %>% 
  mutate(
    param = factor(Parameter, 
                   levels = c("b0", "b1", "b22", "b23", "b3",
                              "sigma.e", "sigma.ui", "sigma.us",
                              "b1s", "b3s", "gamma0", "gamma1"),
                   labels = c(ex(beta[0]), ex(beta[1]), ex(beta[22]), ex(beta[23]), ex(beta[3]),
                              ex(sigma[epsilon]), ex(sigma[U*i]), ex(sigma[U*s]),
                              ex(beta[1*"S"]), ex(beta[3*"S"]), ex(gamma[0]), ex(gamma[1]))
    )
  )

xints <- all_results2 %>% distinct(param)
xints$xint <- c(40, -10, 5, 15, 0.1, 1.5, 3, 2, -0.3, 0.05, 1, 1)

all_results2 %>% 
  ggplot(aes(x = Estimate, fill = id)) + 
  geom_density(alpha = .4) + 
  geom_vline(data = xints, aes(xintercept = xint), alpha = .8, lty = 3, colour = "black") + 
  facet_wrap(~param, scales = "free", labeller = label_parsed, ncol = 2, nrow = 6) + 
  labs(fill = "") + 
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black")
  )
ggsave("./JM-sims-plots/SampleSizeRE.png")

save(all_results2, file = "./JMRE.RData")
save(simstudy, file = "./JMSimStudy.RData")





  