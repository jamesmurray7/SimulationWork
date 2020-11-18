#' ############
#' Expanding Simulating_LMM.R 
#' 
#' This is done by inclusion of extra covariates
#' and more simulated covariate processes
#' (N = 3000, m_i = 6)
#' ############

dev.off()
rm(list=ls())
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(broom.mixed)

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
sigma = 3 # Random effects
epsilon = 1.5 # Measurement error

# Single-run ----

nsubj = N/m
id <- rep(1:nsubj, each = m)
x1 <- rep(rbinom(nsubj, size = 1, prob = .5), each = m)
x2 <- gl(3, m, N)
x3 <- floor(rep(rnorm(nsubj, 60, 10), each = m))
re <- rep(rnorm(nsubj, 0, sigma), each = m)
eps <- rnorm(N, 0, epsilon)

# Formulate y
y = beta_0 + beta_1 * x1 + beta_21 * (x2 == 2) + beta_22 * (x2 == 3) + 
  beta_3 *x3 + re + eps

dat <- data.frame(id, x1, x2, x3, re, eps, y)

X <- model.matrix(~x1 + x2 + x3)
B <- matrix(c(beta_0, beta_1, beta_21, beta_22, beta_3), nrow = 1)

Y = X %*% t(B) + re + eps

fit <- lmer(y ~ x1 + x2 + x3 + (1|id), data = dat)
summary(fit)

# Functionise ~
lmm_sim = function(N = 3000, m = 6,
                   beta_0 = 40, beta_1 = -10,
                   beta_21 = 5, beta_22 = 15, beta_3 = 0.15,
                   sigma = 3, epsilon = 1.5){
  
  # Check we get a whole number of subjects...
  if(str_detect(as.character(N/m),"\\.")){
    stop("Please use divisible N / m")
  }
  
  nsubj = N/m
  id <- rep(1:nsubj, each = m)
  
  x1 <- rep(rbinom(nsubj, size = 1, prob = .5), each = m) # Treatment
  x2 <- gl(3, m, N) # 
  x3 <- floor(rep(rnorm(nsubj, 60, 10), each = m))
  re <- rep(rnorm(nsubj, 0, sigma), each = m)
  eps <- rnorm(N, 0, epsilon)
  
  X <- model.matrix(~x1 + x2 + x3)
  B <- matrix(c(beta_0, beta_1, beta_21, beta_22, beta_3),
              nrow = 1)
  
  Y <- X %*% t(B) + re + eps
  
  return(
    data.frame(
      id, x1, x2, x3, re, eps, Y
    )
  )
}

dat = replicate(1000, lmm_sim(), simplify = F)

lmm_model_fit <- function(x){
  pb$tick()$print()
  fit <- lmer(Y ~ x1 + x2 + x3 + (1|id), data = x)
}

models <- map(dat, lmm_model_fit)

models %>% 
  map_df(~ tidy(.x, effects = "fixed")) %>% 
  filter(str_detect(term, "^x")) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(colour = "grey20", alpha = .2) + 
  facet_wrap(~term, scales = "free")

# Quick investigation into sigma and eps ----------------------------------

testing <- crossing(
  sigma = c(1, 5, 10),
  eps = c(1, 5, 10)
) %>% 
  mutate(id = glue::glue("sigma = {sigma}, epsilon = {eps}"),
         r = row_number()) %>% 
  group_by(id) %>% 
  nest() %>%
  rename(df=data)

empty <- rep(NA, 500)

lmm_wrap <- function(df){
  sig <- df$sigma
  eps <- df$eps
  
  message("\n---\n", round(df$r/nrow(testing) * 100, 1), "%\n")
  
  beta <- data.frame(b0=empty, b1 = empty, b21 = empty, b22 = empty, b3 = empty)
  for(i in 1:500){
    this.dat <- lmm_sim(sigma = sig, epsilon = eps)
    fit <- lmer(Y ~ x1 + x2 + x3 + (1|id), data = this.dat)
    beta[i,] <- getME(fit, "beta")
    
    # Print progress
    if(i %/% 100 > 0 && i %% 100 == 0){
      message("\n---\nIteration ", i, " done")
    }
  }
  return(beta)
}

testing_results <- testing %>% 
  mutate(betas = map(df, lmm_wrap))

testing_results2 <- testing_results %>% 
  ungroup %>% 
  unnest(c(df, betas)) %>% 
  arrange(sigma, eps) %>% 
  mutate(id = fct_inorder(id)) %>% 
  select(-sigma, -eps, -r) %>% 
  gather("beta", "estimate", -id) 

testing_results2 %>% 
  ggplot(aes(x = estimate, colour = id)) + 
  geom_density() + 
  facet_wrap(~beta, scales = "free") + 
  labs(colour = "",
       x = "Estimate", y = "Density",
       title = expression("Effect of differing magnitudes "*sigma*" and "*epsilon),
       subtitle = "500 simulations")
ggsave("./LMM_Epsilon_Sigma.png")


