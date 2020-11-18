rm(list=ls())
library(tidyverse)
library(survival)

# Simulating as per CRAN.R-Project.org link at foot -----------------------

exp_bh <- function(t, lambda = 0.5) lambda * t ^ 0  # As constant
weibull_bh <- function(t, lambda = 0.5, nu = 1.5) lambda * nu * t ^ (nu - 1)

curve(exp_bh, from = 0, to = 5, lty = 1, ylim = c(0, 2), 
      ylab = expression(h[0](t)), xlab = "t")
curve(weibull_bh, from = 0, to = 5, lty = 1, add = T)
legend("topleft", lty = 1:2, legend = c("Exp", "Weibull"),
       title = "Baseline hazard", bty = "n")

sim_data <- function(dataset, n, baseline, params = list(), coveff = -0.50){
  # Simulating a treatment effect
  x <- rbinom(n , 1, .5)
  # Draw U ~ Unif(0,1)
  U <- runif(n)
  # Simulate survival times depending on the BL hazard (Exp/Weibull)
  if(baseline == "Exp"){
    t <- -log(U)/(params$lambda * exp(coveff * x))
  }else if(baseline == "Weibull"){
    t <- (-log(U)/(params$lambda * exp(coveff * x))) ^ (1/params$nu)
  }else{
    stop("Baseline not one of 'Exp' 'Weibull', stopping...")
  }
  t <- abs(t)
  
  # Make event indicator
  d <- as.numeric(t < 5)
  t <- pmin(t, 5)
  # Return data frame
  data.frame(
    dataset = dataset,
    x = x,
    t = t,
    d = d,
    n = n,
    baseline = baseline,
    stringsAsFactors = F
  )
}

sim_data(dataset = 1:100, n = 50, baseline = "Exp", params = list(lambda = .5), 
         coveff = -0.5)

data <- list()
data[["n=50, baseline = Exp"]] <- lapply(
  X = 1:100,
  FUN = sim_data,
  n = 50,
  baseline = "Exp", params = list(lambda = 0.5),
  coveff = -0.5
)
data[["n=250, baseline = Exp"]] <- lapply(
  X = 1:100,
  FUN = sim_data,
  n = 250,
  baseline = "Exp", params = list(lambda = 0.5),
  coveff = -0.5
)
data[["n=50, baseline = Weibull"]] <- lapply(
  X = 1:100,
  FUN = sim_data,
  n = 50,
  baseline = "Weibull", params = list(lambda = 0.5, nu = 1.5),
  coveff = -0.5
)
data[["n=250, baseline = Weibull"]] <- lapply(
  X = 1:100,
  FUN = sim_data,
  n = 250,
  baseline = "Weibull", params = list(lambda = 0.5, nu = 1.5),
  coveff = -0.5
)

cox_res <- list()

get_cox <- function(idx){
  fits <- c()
  for(i in 1:length(data[[idx]])){
    dat <- data[[idx]][[i]]
    fit <- coxph(Surv(dat$t, dat$d) ~ dat$x)
    out <- round(as.numeric(fit$coefficients), 4)
    fits <- c(fits, out)
  }
  return(fits)
}

cox_res <- data.frame(n = c(50, 50, 250, 250),
                      blhaz = rep(c("Exp", "Weibull"), 2),
                      idx = 1:4)

cox_res2 <- cox_res %>% 
  mutate(
    lab = paste0("n = ", n, ", baseline hazard: ", blhaz),
    betas = map(idx, get_cox)) %>% 
  group_by(lab) %>% 
  unnest(betas) %>% 
  ungroup
  
cox_res2 %>% 
  ggplot(aes(x = betas)) + 
  geom_density() + 
  facet_wrap(~lab)

# BL Hazard of Weibull appears to perform well from this 


# Simulating Weibull only -------------------------------------------------

# Baseline hazard: Weibull

sim_Weib <- function(N, lambda, nu, beta, rateC){
  # N Bernoulli trials (below is same as rbinom)
  x <- sample(x = c(0,1), size = N, replace = T, prob = c(.5, .5))
  
  # Weibull event times
  U <- runif(N)
  tt <- (-log(U)/(lambda * exp(beta * x))) ^ (1/nu)
  
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
  dat <- sim_Weib(1000, 0.02, 1, -0.5, 0.001)
  fit <- coxph(Surv(time, status) ~ x, data = dat)
  beta.hat[k] <- fit$coef
}

mean(beta.hat)
sd(beta.hat)/sqrt(length(beta.hat))


# Function to read-in off data-set
test_set <- crossing(
  N = c(50, 250, 500, 1000), # Num simulated data points
  lambda = c(0.01, 0.1, 0.25, 0.5), # Rate of event
  nu = 1,
  beta = -0.5,
  rateC = c(0.001, 0.01, 0.1, 0.25), # Rate of censoring
  num_iter = c(100, 500, 1000) # Number of simulations
)

test_set2 <- test_set %>% 
  mutate(
    id = glue::glue("N: {N}, lambda: {lambda}, censor rate: {rateC}, {num_iter} iterations")
  ) %>% 
  group_by(id) %>% 
  nest() %>% 
  rename(df=data) %>% 
  ungroup()

sim_weib_wrap <- function(df){
  N = df$N # Number of simulated failure times
  lambda = df$lambda # Rate of failure times
  nu = df$nu # Shape of Weibull function
  beta = df$beta # Target coefficient
  rateC = df$rateC # Rate of censoring
  num_iter = df$num_iter # Number of simulations
  
  cat(N, lambda, nu, beta, rateC,"\n")
  
  betas <- c()
  for(i in 1:num_iter){
    dat <- sim_Weib(N = N,
                    lambda = lambda,
                    nu = nu,
                    beta = beta,
                    rateC = rateC)
    fit <- coxph(Surv(time, status) ~ x, data = dat)
    betas <- c(betas, fit$coef)
    
    # Print progress
    if(i %/% 100 > 0 && i %% 100 == 0){
      message("\n---\nIteration ", i, " done")
    }
  }
  return(betas)
}

test_set2$betas <- map(test_set2$df, sim_weib_wrap)

test_results <- test_set2 %>% 
  group_by(id) %>% 
  unnest(cols = c(df, betas)) %>% 
  ungroup


# Plots -------------------------------------------------------------------

# Effect of increasing N (number sim. data points) ----
test_results %>% 
  filter(
    lambda == 0.01,
    rateC == 0.1
  ) %>% 
  arrange(N, num_iter) %>% 
  mutate(plot_id = glue::glue("Number of simulated event-times: {N},\nnumber of simulations: {num_iter}"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = betas)) + 
  geom_density(colour = "grey20", alpha = .25) + 
  geom_vline(xintercept = -0.5, lty = 5, alpha = .5, colour = "blue") +
  facet_wrap(~plot_id, scales = "free", ncol = 3) + 
  labs(x = expression(hat(beta)), title = "Weibull",
       subtitle = "lambda = 0.01, censor rate = 0.1") + 
  theme_light()
# Save
ggsave("H:/R work/Weibull simulation plots/N_num_iter.png")

# Effect of increasing failure rate ----
test_results %>% 
  filter(
    N == 1000,
    num_iter == 1000,
    rateC == 0.1
  ) %>% 
  arrange(lambda) %>% 
  mutate(plot_id = glue::glue("Lambda: {lambda}"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = betas)) + 
  geom_density(colour = "grey20", alpha = .25) + 
  geom_vline(xintercept = -0.5, lty = 5, alpha = .5, colour = "blue") +
  facet_wrap(~plot_id, scales = "free", ncol = 3) + 
  labs(x = expression(hat(beta)), title = "Weibull",
       subtitle = "Number simulated times = 1000, num simulations = 1000, censor rate = 0.1") + 
  theme_light()
ggsave("H:/R work/Weibull simulation plots/lambda.png")

# Effect of increasing censor rate ----
test_results %>% 
  filter(
    N == 1000,
    num_iter == 1000,
  ) %>% 
  arrange(lambda, rateC) %>% 
  mutate(plot_id = glue::glue("Lambda: {lambda}, C = {rateC}"),
         plot_id = fct_inorder(plot_id)) %>% 
  ggplot(aes(x = betas)) + 
  geom_density(colour = "grey20", alpha = .25) + 
  geom_vline(xintercept = -0.5, lty = 5, alpha = .5, colour = "blue") +
  facet_wrap(~plot_id, scales = "free", ncol = 3) + 
  labs(x = expression(hat(beta)), title = "Weibull",
       subtitle = "Number simulated times = 1000, num simulations = 1000") + 
  theme_light()

ggsave("H:/R work/Weibull simulation plots/lambda_rateC.png")

# https://cran.r-project.org/web/packages/rsimsum/vignettes/B-relhaz.html
# https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring

