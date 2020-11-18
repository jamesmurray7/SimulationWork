#' ############
#' Expanding Simulating_LMM_p2.R 
#' 
#' Covariates kept the same, 
#' Now introducing MCAR missingness
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

# Scenario we wish to simulate:
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

# Cast to model matrix
X <- model.matrix(~x1 + x2 + x3)
B <- matrix(c(beta_0, beta_1, beta_21, beta_22, beta_3), nrow = 1)

Y = X %*% t(B) + re + eps


dat <- data.frame(id, x1, x2, x3, re, eps, Y, mmm)
dat$Y <- ifelse(dat$mmm < 0.33, NA, dat$Y) # 33% Missingness

fit <- lmer(Y ~ x1 + x2 + x3 + (1|id), data = dat)
summary(fit)

# Functionise ~
lmm_sim_MCAR = function(N = 3000, m = 6, MCAR,
                   beta_0 = 40, beta_1 = -10,
                   beta_21 = 5, beta_22 = 15, beta_3 = 0.15,
                   sigma = 2.5, epsilon = 2.5){
  
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
  
  mm <- runif(N)
  out <- data.frame(id, x1, x2, x3, re, eps, Y, mm)
  
  out$Y <- ifelse(out$mm < 0.33, NA, out$Y)
  
  return(out)
}

dat = replicate(1000, lmm_sim_MCAR(MCAR = 0.33), simplify = F)

lmm_model_fit <- function(x){
  pb$tick()$print()
  fit <- lmer(Y ~ x1 + x2 + x3 + (1|id), data = x)
}
pb <- progress_estimated(1000)
models <- map(dat, lmm_model_fit)

models %>% 
  map_df(~ tidy(.x, effects = "fixed")) %>% 
  filter(str_detect(term, "^x")) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(colour = "grey20", alpha = .2) + 
  facet_wrap(~term, scales = "free")

# Quick investigation into degree of missingness and errors ---------------

testing <- crossing(
  MCAR = c(0.1, 0.2, 0.5, 0.75),
  sigma = c(1, 5, 10),
  eps = c(1, 5, 10)
) %>% 
  mutate(id = glue::glue("MCAR = {MCAR}, sigma = {sigma}, epsilon = {eps}"),
         r = row_number()) %>% 
  group_by(id) %>% 
  nest() %>%
  rename(df=data)

empty <- rep(NA, 500)

lmm_wrap <- function(df){
  mcar <- df$MCAR
  sig <- df$sigma
  eps <- df$eps
  
  message("\n---\n", round(df$r/nrow(testing) * 100, 1), "%\n")
  
  beta <- data.frame(b0=empty, b1 = empty, b21 = empty, b22 = empty, b3 = empty)
  
  pb <- progress::progress_bar$new(total = 500)
  
  for(i in 1:500){
    this.dat <- lmm_sim_MCAR(MCAR = mcar, sigma = sig, epsilon = eps)
    suppressMessages(fit <- lmer(Y ~ x1 + x2 + x3 + (1|id), data = this.dat))
    beta[i,] <- getME(fit, "beta")
    
    # Print progress
    pb$tick()
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
  arrange(MCAR, sigma, eps) %>% 
  mutate(id = fct_inorder(id)) %>% 
  select(-sigma, -eps, -r, -MCAR) %>% 
  gather("beta", "estimate", -id) 

mcar_plot <- function(mcar){
  string <- paste0("MCAR = ", mcar)
  testing_results2 %>% 
    filter(str_detect(id, string),
           beta != "b0") %>% 
    mutate(id = str_extract(id, "sigma.*$")) %>% 
    ggplot(aes(x = estimate, colour = id)) + 
    geom_density() + 
    facet_wrap(~beta, scales = "free") + 
    labs(colour = "", x = "", y = "",
         title = string)
  
}

p1 <- mcar_plot(0.1)
p2 <- mcar_plot(0.2)
p3 <- mcar_plot(0.5)
p4 <- mcar_plot(0.75)

ggpubr::ggarrange(p1,p2,p3,p4,
                  common.legend = T, legend = "bottom",
                  ncol = 2, nrow = 2)
ggsave("./MCAR_LMM.png")

