# For Bayesian modeling
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stroop_m1 <- stan_model("model1.stan")

fit_m1 <- sampling(
  stroop_m1,
  data = stan_dat,
  iter = 2000,
  warmup = 500,
  chains = 3,
  cores = 3,
  seed = 2
)

pairs(fit_m1, pars = c("mu_beta_con", "sigma_con"))