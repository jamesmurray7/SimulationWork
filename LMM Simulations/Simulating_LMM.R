rm(list=ls())
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(broom.mixed)

# Random slope ----

nstand = 5
nplot = 4
mu = 10 # TRUE mean
sds = 2 # shared random effects
sd = 1 # epsilon error

stand = rep(LETTERS[1:nstand], each = nplot)
plots = letters[1:(nstand*nplot)]

set.seed(16)
# Random effects
raneff =  rnorm(nstand, 0, sds)
raneff = rep(raneff, each = nplot)
# Epsilon (subject-level) error
eps = rnorm(nstand * nplot, 0, sd)

dat = data.frame(stand, raneff, plots, eps)
# Fit data y
dat$y = with(dat, mu + raneff + eps)
# LMER model
fit = lmer(y ~ 1 + (1|stand), data = dat)
fit # 10.57

# Functionise
lmm_sim = function(ngroup = 5,
                   nsubj = 4,
                   mu = 10,
                   sigma = 2,
                   sd = 1){
  raneff = rep(rnorm(ngroup, 0 , sigma), each  = nsubj)
  group = rep(LETTERS[1:ngroup], each = nsubj)
  eps = rnorm(nstand*nplot, 0, sd)
  resp = mu + raneff + eps
  dat = data.frame(group, resp)
  lmer(resp ~ 1 + (1|group), data = dat)
}
lmm_sim()

sims = replicate(100, lmm_sim(), simplify = F)

sims %>% 
  map_df(~tidy(.x, effects = "fixed")) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .15)

# Actual fixed-effects ----
rm(list = ls())

ngroup = 5 # Five groups of subjects
nsubj = 4  # Four subjects per group
b0 = -1   # Mean response
b1 = .005 # Change in mean response for 1 unit change in x1
b2 = .100 # Change in mean response for 1 unit change in x2
sds = 2 # SD of random effects
sd = 1 # subject-level SD (epsilon error)

set.seed(16)
group = rep(LETTERS[1:ngroup], each = nsubj)
raneff = rep(rnorm(ngroup, 0, sds), each = nsubj)
eps = rnorm(ngroup * nsubj)

x1 = rep(runif(ngroup, 1000, 1500), each = nsubj) # Across group
x2 = runif(ngroup * nsubj, 2, 75) # Subject-specific

y = b0 + b1 * x1 + b2 * x2 + raneff + eps
lmer(y ~ x1 + x2 + (1|group))

# Functionise ~
lmm_sim = function(ngroup = 5, nsubj = 4, b0 = -1, b1 = 0.005, b2 = .1,
                   sds = 2, sd = 1){
  group = rep(LETTERS[1:ngroup], each = nsubj)
  raneff = rep(rnorm(ngroup, 0, sds), each = nsubj)
  eps = rnorm(ngroup * nsubj, 0, sd = sd)
  
  x1 = rep(runif(ngroup, 1000, 1500), each = nsubj)
  x2 = runif(ngroup * nsubj, 2, 75)
  
  y = b0 + b1 * x1 + b2 * x2 + raneff + eps
  dat = data.frame(group, x1, x2, y)
  lmer(y ~ x1 + x2 + (1|group))
}

sim = replicate(1000, lmm_sim(), simplify = F)
sim %>% 
  map_df(~ tidy(.x, effects = "fixed")) %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(colour = "grey20", alpha = .2) + 
  facet_wrap(~term, scales = "free")
