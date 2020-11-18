rm(list=ls())
library(tidyverse)
theme_set(theme_light())
library(broom)
set.seed(16)

# Example 1 ----

#' Want to compare two groups, 
#' group one mean 5, group two is 2 lower
#' sd of residuals 2 

ngroup = 2
nrep = 10
b0 = 5
b1 = -2
sd = 2

group = rep(c("group1", "group2"), each = nrep)

# Simulate the random errors
eps = rnorm(n = nrep * ngroup, mean = 0, sd = sd)

# Generate Y values
y = b0 + b1 * (group == "group2") + eps

df = data.frame(group, y)
# Fit one model
fit = lm(y ~ group, df)
summary(fit)

# Functionise the simulated data ----
simfn = function(nrep = 10, b0 = 5, b1 = -2, sigma = 2){
  ngroup = 2
  groups = rep(c("group1", "group2"), each = nrep)
  eps = rnorm(n = ngroup * nrep, mean = 0, sd = sigma)
  y = b0 + b1 * (group == "group2") + eps
  simdat = data.frame(group, y)
  simfit = lm(y ~ group, simdat)
  simfit
}
simfn()
simfn(sigma = .1)

# Repeat this simulation many times
sims = replicate(n = 1000, simfn(), simplify = F)
sims[[5]]

sims %>% 
  map_df(tidy) %>% 
  filter(term == "groupgroup2") %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .4) + 
  geom_vline(xintercept = -2)

# Plotting standard deviation of residuals
sims %>% 
  map_dbl(~summary(.x)$sigma) %>% 
  data.frame(sigma = .) %>% 
  ggplot(aes(x = sigma)) + 
  geom_density(fill = "grey20", alpha = .4) + 
  geom_vline(xintercept = 2) + 
  labs(x = expression(sigma))

# Is effect significant?
sims %>% 
  map_df(tidy) %>% 
  filter(term == "groupgroup2") %>%
  summarise(med = median(p.value), mean = mean(mean(p.value)))


# Example 2 ----
rm(list=ls())
# Group 1 mean 10, group 2 mean + 3, -.25 per unit age

ngroup = 2
nrep = 10
b0 = 10
b1 = 3
b2 = -0.25
sd = 2.5

# Covariates
group = rep(c("g1", "g2"), each = nrep)
age = floor(rnorm(ngroup*nrep, 50, 5))
# Error
eps = rnorm(ngroup*nrep, 0, sd)
# Y-values
y = b0 + b1 * (group == "g2") + b2 * age + eps

df = data.frame(
  group = group,
  age = age,
  y = y
)

fit = lm(y ~ age + group, df)
summary(fit)

# Functionise and wrap
sim = function(nrep = 10, ngroup = 2,
               b0 = 10, b1 = 3, b2 = -.25, 
               sigma = 2.5){
  group = rep(c("g1", "g2"), each = nrep)
  age = floor(rnorm(ngroup*nrep, 50, 5))
  eps = rnorm(ngroup*nrep, 0, sigma)
  y = b0 + b1 * (group == "g2") + b2 * age + eps
  simdf = data.frame(group = group, age = age, y = y)
  simfit = lm(y ~ age + group, simdf)
  simfit
}

sims = replicate(1000, sim(), simplify = F)
sims_lowsigma = replicate(1000, sim(sigma = 1), simplify = F)  
sims_morerep  = replicate(1000, sim(nrep = 100), simplify = F)  
  
sims_morerep %>% 
  map_df(tidy) %>% 
  filter(str_detect(term, "g2$|age")) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .25) + 
  facet_wrap(~term, scales = "free")

sims_morerep %>% 
  map_dbl(~summary(.x)$sigma) %>% 
  tibble() %>% 
  ggplot(aes(x = .)) + 
  geom_density(fill = "grey20")

  
  
  