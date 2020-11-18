#' ######################################
#' Builds on work in simulating survival times .R
#' by introducing Gompertz simulation and repeating process
#' (essentially)
#' #####################################


# Prerequisites -----------------------------------------------------------
rm(list=ls())
library(tidyverse)
theme_set(theme_light())
library(survival)

# Plotting out what baseline hazards look like ----------------------------

exp_bh <- function(t, lambda = 0.5) lambda * t ^ 0  # As constant
weibull_bh <- function(t, lambda = 0.5, nu = 1.5) lambda * nu * t ^ (nu - 1)
gompertz_bh <- function(t, lambda = 0.5, alpha = .2) lambda * exp(alpha * t)

curve(exp_bh, from = 0, to = 5, lty = 1, ylim = c(0, 2), 
      ylab = expression(h[0](t)), xlab = "t")
curve(weibull_bh, from = 0, to = 5, lty = 2, add = T)
curve(gompertz_bh, from = 0, to = 5, lty = 3, add = T)
legend("topleft", lty = 1:3, legend = c("Exp", "Weibull", "+Gompertz"),
       title = expression("baseline hazard "*lambda*" = 0.5"), bty = "n")

# Let's look at how changing shape and scale affects Gompertz for t [0,20]

gz_haz <- function(lambda, alpha){
  # Lambda - scale (rate); alpha - shape
  t <- 0:20
  g <- lambda * exp(alpha * t)
  plot(t,g, type = "l")
  text(x = 5, y = max(g), paste0("shape = ", alpha, " scale = ", lambda))
  return(g)
}

test_set <- crossing(
    lambda = c(0.01, 0.1, 0.25, 0.5),
    alpha = c(0.01, 0.1, 0.25, 0.5)
  ) %>% 
  mutate(id = glue::glue("Lambda = {lambda}, alpha = {alpha}"))

test_results <- test_set %>% 
  mutate(
    gompertz = map2(lambda, alpha, gz_haz)
  )

test_results <- test_results %>% 
  group_by(id) %>% 
  unnest(gompertz) %>% 
  mutate(t = 0:20) %>%
  ungroup

test_results %>% 
  ggplot(aes(x = t, y = gompertz)) + 
  geom_line(lwd = 1.5, colour = "magenta") + 
  facet_wrap(~id, scales = "free", nrow = 4, ncol = 4) + 
  theme_light()

# Simulating Gompertz times -----------------------------------------------

# Gompertz Baseline hazard - one binomial covariate.

sim_Gomp <- function(N, lambda, alpha, beta, rateC){
  # N Bernoulli trials (below is same as rbinom)
  x <- sample(x = c(0,1), size = N, replace = T, prob = c(.5, .5))
  
  # Gompertz event times
  U <- runif(N)
  tt <- 1/alpha * log(1 - (alpha * log(U))/(lambda * exp(beta * x)))
  
  # Censoring times
  C <- rexp(N, rate = rateC)
  
  # Follow-up times and event indicators
  time <- pmin(tt, C)
  status <- as.numeric(tt <= C)
  
  # Return data set
  data.frame(
    id = 1:N,
    x = x,
    time = time,
    status = status
  )
}

beta.hat <- c()
for(k in 1:1e3){
  dat <- sim_Gomp(100, 0.01, 0.25, -2, 0.01)
  fit <- coxph(Surv(time, status) ~ x, data = dat)
  beta.hat[k] <- fit$coef
  if(k %/% 100 == 10) message(k,"\n")
}

mean(beta.hat) # mean beta.hat
sd(beta.hat)/sqrt(length(beta.hat)) # SE[beta.hat]

# Estimation of three covariates ------------------------------------------
coefficients(lm(runif(300) ~ gl(3,100)))

simmv_Gomp <- function(N, lambda, alpha, beta1, beta2, beta3, beta4, rateC){

  # Beta-vector
  # beta1: Treatment, beta23: x2-2 x2-3; beta4: Age
  betas <- matrix(c(beta1, 0, beta2, beta3, beta4), ncol = 1)
  
  x1 <- rbinom(N, 1, 0.5)
  x2 <- gl(3, 1, length = N)
  x3 <- rnorm(N, 70, 10)
  
  dat <- data.frame(x1,x2,x3)
  dat_mm <- model.matrix(~ . - 1, dat)
  
  # Gompertz event times
  U <- runif(N)
  tt <- 1/alpha * log(1 - (alpha * log(U))/(lambda * exp(dat_mm %*% betas)))
  
  # Censoring?
  C <- rexp(N, rateC)
  
  # Follow-up times and event indicators
  time <- pmin(tt, C)
  status <- as.numeric(tt <= C)
  
  # Return data set
  data.frame(
    id = 1:N,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    time = time,
    status = status
  )
  
}

head(simmv_Gomp(N = 100,
                lambda = 0.01,
                alpha = 0.25,
                beta1 = -2,
                beta2 = 0.5,
                beta3 = 0.75,
                beta4 = 0.05,
                rateC = 0.001))

dat <- simmv_Gomp(N = 100,
                  lambda = 0.01,
                  alpha = 0.25,
                  beta1 = -2,
                  beta2 = 0.5,
                  beta3 = 0.75,
                  beta4 = 0.05,
                  rateC = 0.001)

fit <- coxph(Surv(time, status) ~ x1 + x2 + x3, data = dat)
coef(fit)

beta1 <- beta2 <- beta3 <- beta4 <- c()

for(i in 1:1000){
  dat <- simmv_Gomp(N = 100,
                   lambda = 0.01,
                   alpha = 0.25,
                   beta1 = -2,
                   beta2 = 0.5,
                   beta3 = 0.75,
                   beta4 = 0.05,
                   rateC = 0.001)
  fit <- coxph(Surv(time, status) ~ x1 + x2 + x3, data = dat)
  beta1 <- c(beta1, coef(fit)[1])
  beta2 <- c(beta2, coef(fit)[2])
  beta3 <- c(beta3, coef(fit)[3])
  beta4 <- c(beta4, coef(fit)[4])
  # Print progress
  if(i %/% 100 > 0 && i %% 100 == 0){
    message("\n---\nIteration ", i, " done")
  }
}

mnse <- function(x){
  mn <- mean(x)
  se <- sd(x)/sqrt(length(x))
  return(paste0(mn, " (", se, ")"))
}

betas <- list(beta1, beta2, beta3, beta4)
lapply(betas, mnse)


# Seeing how changing N, shape, scale of Gompertz changes affects  --------

# Function to read-in off data-set

head(simmv_Gomp(N = 100,
                  lambda = 0.01,
                  alpha = 0.25,
                  beta1 = -2,
                  beta2 = 0.5,
                  beta3 = 0.75,
                  beta4 = 0.05,
                  rateC = 0.001))


test_set <- crossing(
  N = c(50, 250, 500), # Num simulated data points
  lambda = c(0.01, 0.1), # Rate of event
  alpha = c(0.1, 0.25), # Shape of Gompertz
  beta1 = -2, # Log-odds associated with binary covariate present.
  beta2 = c(0.5, 1.5), # Log-odds associated with 2nd level of factor covariate
  beta3 = c(0.75, 2.5), # Log-odds associated with 3rd level of factor covariate
  beta4 = 0.05, # Log-odds per unit increase of continuous variable
  rateC = c(0.001, 0.01, 0.25), # Rate of censoring
  num_iter = c(100, 500, 1000) # Number of simulations
)

test_set2 <- test_set %>% 
  mutate(
    id = glue::glue(
      "N: {N}, lambda: {lambda}, alpha = {alpha},
      beta2 = {beta2}, beta3 = {beta3}, censor rate: {rateC}, {num_iter} iterations"
    )
  ) %>% 
  mutate(r = row_number()) %>% 
  group_by(id) %>%  
  nest() %>% 
  rename(df=data) %>% 
  ungroup() 
  

Gompertz_wrap <- function(df){
  N = df$N # Number of simulated failure times
  lambda = df$lambda # Rate of failure times
  alpha = df$alpha # Shape of Weibull function
  rateC = df$rateC # Rate of censoring
  num_iter = df$num_iter # Number of simulations
  
  # Target coefficients
  beta1 = df$beta1
  beta2 = df$beta2
  beta3 = df$beta3
  beta4 = df$beta4
  
  message(round(df$r/49 * 100, 1), "%----\n")
  
  cat("\nN =", N,
      " lambda =", lambda,
      " alpha =", alpha, 
      " beta2 =", beta2,
      " beta3 =", beta3,
      " censor rate = ", rateC,
      "beginning ", num_iter, "iterations...\n")
  
  empty <- rep(NA, num_iter)
  betas <- data.frame(beta1 = empty, beta2 = empty, beta3 = empty, beta4 = empty)
  for(i in 1:num_iter){
    dat <- simmv_Gomp(N = N,
                      lambda = lambda,
                      alpha = alpha,
                      beta1 = beta1,
                      beta2 = beta2,
                      beta3 = beta3,
                      beta4 = beta4,
                      rateC = rateC)
    fit <- coxph(Surv(time, status) ~ x1 + x2 + x3, data = dat)
    betas$beta1[i] <- coef(fit)[1]
    betas$beta2[i] <- coef(fit)[2]
    betas$beta3[i] <- coef(fit)[3]
    betas$beta4[i] <- coef(fit)[4]
    # Print progress
    if(i %/% 100 > 0 && i %% 100 == 0){
      message("\n---\nIteration ", i, " done")
    }
  }
  return(betas)
}

test_set2$betas <- map(test_set2$df, Gompertz_wrap)

test_results <- test_set2 %>% 
  group_by(id) %>% 
  unnest(cols = df) %>% 
  rename_at(vars(beta1:beta4), ~ paste0("true", .x)) %>% 
  unnest(cols = c(betas)) %>% 
  ungroup


# Plots -------------------------------------------------------------------

# Effect of increasing N (number sim. data points) ----

# Beta_1 (Binary)

test_results %>% 
  filter(truebeta2 == 0.5, truebeta3 == 0.75,
         lambda == 0.01, alpha == 0.1, rateC == 0.001) %>% 
  arrange(N, num_iter) %>% 
  mutate(plot_id = glue::glue("{N} survival times\n{num_iter} simulations"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = beta1)) + 
  geom_vline(xintercept = -2, colour = "red", lty =3) + 
  geom_density(fill = "grey20", alpha = .2) + 
  facet_wrap(~plot_id) + 
  labs(x = expression(beta[1]),
       y = "Density",
       title = "Binary covariate",
       subtitle = expression(lambda*"=0.01; "*alpha*"=0.1; censor rate=0.001"))
ggsave("./Gompertz simulation plots/N_num_iter.png")

beta23 <- test_results %>% 
  filter(N == 250, num_iter == 500,
         lambda == 0.01, alpha == 0.1, rateC == 0.001) %>% 
  arrange(truebeta2, truebeta3) %>% 
  mutate(plot_id = glue::glue("beta_2 = {truebeta2}, beta_3 = {truebeta3}"),
         plot_id = fct_inorder(plot_id))

beta23 %>% 
  select(plot_id, beta2, beta3) %>% 
  gather("beta", value = "estimate", -plot_id) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .2) + 
  facet_wrap(plot_id ~ beta) +
  labs(x = expression(beta),
       y = "Density",
       title = "Three-levelled categorical covariate",
       subtitle = expression(lambda*"=0.01; "*alpha*"=0.1; censor rate=0.001,
                             500 simulations of 250 survival times"))
ggsave("./Gompertz simulation plots/Magnitude_categorical_covariates.png")  
  

# Effect of varying shape and scale parameters in Gompertz ----

test_results %>% 
  filter(N == 250, num_iter == 500,
         truebeta2 == 1.5, truebeta3 == 2.5,
         rateC == 0.001) %>% 
  arrange(lambda, alpha) %>% 
  mutate(plot_id = glue::glue("Shape: {lambda}; scale: {alpha}"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = beta1)) + 
  geom_vline(xintercept = -2, colour = "red", lty = 3) + 
  geom_density(fill = "grey20", alpha = .2) + 
  facet_wrap(~plot_id)

# Effect of increasing censor rate ----

test_results %>% 
  filter(truebeta2 == 1.5, truebeta3 == 2.5,
         lambda == 0.01, alpha == .1) %>% 
  arrange(N, num_iter, rateC) %>% 
  mutate(plot_id = glue::glue("{num_iter} simulations\n{N} survival times\nat rate of censoring {rateC}"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = beta1)) + 
  geom_vline(xintercept = -2, colour = "blue", lty = 3, alpha = .5) +
  geom_density() + 
  facet_wrap(~plot_id) +
  labs(x = expression(beta[1])) + 
  theme(strip.text = element_text(size = 7))
  
ggsave("./Gompertz simulation plots/N_num_iter_rateC.png")

# Greater investigation into effect of shape/scale ------------------------

scale_shape <- crossing(
  lambda = c(0.001, 0.01, seq(0.1, 0.5, 0.1)),
  alpha = c(0.01, 0.05, seq(0.1, 0.5, 0.1)),
  rateC = 0.01,
  beta1 = -2, beta2 = 1.5, beta3 = 1.5, beta4 = 0.05,
  N = 250, num_iter = 500
) 

scale_shape2 <- scale_shape %>% 
  mutate(id = glue::glue("lambda: {lambda}, alpha: {alpha}"),
         r = row_number()) %>% 
  group_by(id) %>% 
  nest() %>% 
  rename(df=data) %>% 
  ungroup() 
  
scale_shape2$betas <- map(scale_shape2$df, Gompertz_wrap)  

scale_shape_results <- scale_shape2 %>% 
  unnest(df) %>% 
  rename_at(vars(beta1:beta4), ~ paste0("true", .x)) %>% 
  unnest(betas)

scale_shape_results %>% 
  arrange(lambda, alpha) %>% 
  mutate(plot_id = id,
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = beta1)) + 
  geom_vline(xintercept = -2, colour = "blue", lty = 3) + 
  geom_density() + 
  facet_wrap(~plot_id) + 
  theme(strip.text = element_text(size = 6))
ggsave("./Gompertz simulation plots/scale_shape.png")
