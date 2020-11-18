#' ############
#' Expanding Simulating_LMM_p2 - and _p3.R 
#' 
#' This is done by keeping covariates the same
#' and introducing time covariate (for random slopes).
#' (N = 3000, m_i = 6)
#' ############

# Prerequisites -----------------------------------------------------------
dev.off()
rm(list=ls())
library(MASS)
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(broom.mixed)
library(MASS)

# Fixed-effects and random slope per subject ----

m <- 6 # 6 measurements per subject
N = 3e3 # We want a 'long' data-frame with 3000 measurements.

# Let's think of a scenario we wish to simulate:
# Clinical trial, outcome some diabetes measurement (lower = better)
# Binary: Treatment
# Categorical factor: BMI group at baseline
# Continuous: Age at baseline (mean 60)

# Mean response: 
beta_0 = 40
beta_1 = -10 # Receiving treatment
beta_21 = 5 # Overweight effect
beta_22 = 15 # Obese effect
beta_3 = 0.15 # Effect per unit age
# Error terms
sigma.i = 3 # Random effects - intercept
sigma.s = 2 # Random effects - slope
epsilon = 1.5 # Measurement error

# Single-run ----

nsubj = N/m
id <- rep(1:nsubj, each = m)
x1 <- rep(rbinom(nsubj, size = 1, prob = .5), each = m)
x2 <- gl(3, m, N)
x3 <- floor(rep(rnorm(nsubj, 60, 10), each = m))
time <- rep(0:5, nsubj)
# Measurement error
eps <- rnorm(N, 0, epsilon)
# Random effects
re <-  mvrnorm(nsubj, mu=c(0,0), Sigma=rbind(c(sigma.i, sigma.s),
                                          c(sigma.s, sigma.i)))
re.i <- rep(re[,1], each = m)
re.s <- rep(re[,2], each = m)

# Cast to model matrix
X <- model.matrix(~x1 + x2 + x3)
B <- matrix(c(beta_0, beta_1, beta_21, beta_22, beta_3), nrow = 1)

Y = X %*% t(B) + re.i + re.s * time + eps

dat <- data.frame(id, time, x1, x2, x3, re.i, re.s, eps, Y)

fit <- lmer(Y ~ time + x1 + x2 + x3 + (1+time|id), data = dat, REML = F)
summary(fit)

# Functionise ~
lmm_sim = function(N = 3000, m = 6,
                   beta_0 = 40, beta_1 = -10,
                   beta_21 = 5, beta_22 = 15, beta_3 = 0.15,
                   sigma.i, sigma.s, rho = 0.25,
                   epsilon = 1.5){
  
  # Check we get a whole number of subjects...
  if(str_detect(as.character(N/m),"\\.")){
    stop("Please use divisible N / m")
  }
  
  # Covariance matrix
  
  sigma <- matrix(
    c(sigma.i ^ 2, rho * sigma.i * sigma.s,
      rho * sigma.i * sigma.s, sigma.s ^ 2),
    nrow = 2, byrow = T
  )
  
  stopifnot(eigen(sigma)$values > 0) # positive-definite
  
  nsubj = N/m
  id <- rep(1:nsubj, each = m)
  
  x1 <- rep(rbinom(nsubj, size = 1, prob = .5), each = m) # Treatment
  x2 <- gl(3, m, N) # 
  x3 <- floor(rep(rnorm(nsubj, 60, 10), each = m))
  # Measurement error
  eps <- rnorm(N, 0, epsilon)
  # Random effects
  re <-  mvrnorm(nsubj, mu = c(0, 0), Sigma = sigma)
  re.i <- rep(re[,1], each = m)
  re.s <- rep(re[,2], each = m)
  
  X <- model.matrix(~x1 + x2 + x3)
  B <- matrix(c(beta_0, beta_1, beta_21, beta_22, beta_3),
              nrow = 1)
  
  Y <- X %*% t(B) + re.i + re.s * time + eps
  
  dat <- data.frame(id, time, x1, x2, x3, re.i, re.s, eps, Y)
  
  return(dat)
}

dat = replicate(1000, lmm_sim(sigma.i = 10, sigma.s = 3), simplify = F)

lmm_model_fit <- function(x){
  pb$tick()
  fit <- lmer(Y ~ x1 + x2 + x3 + time + (1|id) + (0 + time|id), data = x,
              control = lmerControl(optimizer = "Nelder_Mead"))
}
pb <- progress::progress_bar$new(total = 1000)

models <- map(dat, lmm_model_fit)

models2 <- models %>% 
  map_df(~ tidy(.x)) %>% 
  filter(str_detect(term, "^x|^sd")) 

models2 %>% 
  mutate(
    term2 = factor(term,
                   levels = c("x1","x22","x23","x3","sd__(Intercept)",
                              "sd__time","sd__Observation"),
                   labels = c(expression(beta[1]),
                              expression(beta[21]),
                              expression(beta[22]),
                              expression(beta[3]),
                              expression(sigma[i]),
                              expression(sigma[t]),
                              expression(epsilon))
    )
  ) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(colour = "grey20", alpha = .2) + 
  facet_wrap(~term2, scales = "free", labeller = label_parsed) + 
  theme(strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_blank())

# Quick investigation into sigma_i and sigma_t and rho --------------------

# check for positive-definiteness
pd_check <- function(sigma.i, sigma.s, rho){
  sigma <- matrix(
    c(sigma.i ^ 2, rho * sigma.i * sigma.s,
      rho * sigma.i * sigma.s, sigma.s ^ 2),
    nrow = 2, byrow = T
  )
  flag <- 0
  if(all(eigen(sigma)$value > 0)) flag <- 1
  return(flag)
}

testing <- crossing(
  sigma.i = c(5, 10, 20),
  sigma.s = c(2, 5, 10),
  rho = c(0.2, 0.6),
  eps = c(2, 5, 10)
) %>% 
  mutate(id = glue::glue("re.i = {sigma.i}, re.s = {sigma.s}, corr = {rho}, epsilon = {eps}"),
         r = row_number(),
         flag = pmap_dbl(list(sigma.i, sigma.s, rho), pd_check)) %>% 
  group_by(id, flag) %>% 
  nest() %>%
  rename(df=data) %>% 
  ungroup


# Begin function
empty <- rep(NA, 250)

lmm_wrap <- function(df){
  sigma.i <- df$sigma.i
  sigma.s <- df$sigma.s
  rho <- df$rho
  eps <- df$eps
  
  # Mainly for debugging
  message(round(df$r/nrow(testing) * 100), "%")
  message("i=",sigma.i," s=", sigma.s, " rho=", rho, " eps=", eps)
  
    # Intialise empty df of parameters
  params <- data.frame(b1 = empty, b21 = empty, b22 = empty, b3 = empty,
                      # SDs
                      sigma.i = empty, sigma.s = empty, eps = empty)
  pb_loop <- progress::progress_bar$new(total = 250)
  for(i in 1:250){
    # Simulate data
    this.dat <- lmm_sim(sigma.i = sigma.i,
                        sigma.s = sigma.s,
                        epsilon = eps, rho = rho)
    # Fit model
    fit <- lmer(Y ~ x1 + x2 + x3 + time + (1|id) + (0 + time|id), data = this.dat,
                control = lmerControl(optimizer = "Nelder_Mead"))
    
    # Get parameter estimates
    this.fit <- tidy(fit)
    this.param <- this.fit %>% 
      filter(str_detect(term, "^x|^sd")) %>% 
      pull(estimate)
    
    params[i, ] <- this.param
    
    # Print progress
    pb_loop$tick()
    if(i %/% 100 > 0 && i %% 100 == 0){
      message("\n---\nIteration ", i, " done")
    }
  }
  return(params)
}

testing_results <- testing %>% 
  mutate(params = map(df, lmm_wrap))

testing_results2 <- testing_results %>% 
  ungroup %>% 
  unnest(df) %>% 
  filter(rho == 0.2, eps %in% c(5, 10)) %>% 
  arrange(sigma.i, sigma.s, rho, eps) %>% 
  select(-sigma.i, -sigma.s, -eps) %>% 
  unnest(params) %>% 
  mutate(id = fct_inorder(id)) %>% 
  select(-flag, -rho, -r) %>% 
  gather("parameter", "estimate", -id) 

testing_results2 %>% 
  mutate(
    parameter2 = factor(parameter,
                   levels = c("b1","b21","b22","b3","sigma.i",
                              "sigma.s","eps"),
                   labels = c(expression(beta[1]),
                              expression(beta[22]),
                              expression(beta[23]),
                              expression(beta[3]),
                              expression(sigma[i]),
                              expression(sigma[t]),
                              expression(epsilon))
    ),
    id = str_remove(id, "corr.*\\, ")
  ) %>% 
  ggplot(aes(x = estimate, colour = id)) + 
  geom_density() + 
  facet_wrap(~parameter2, scales = "free", labeller = label_parsed) + 
  labs(colour = "") + 
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12, colour = "black"),
        strip.background = element_blank())
ggsave("./LMM_Epsilon_Sigma2.png")


