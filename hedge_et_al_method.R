## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

# For easy data manipulation
library(foreach)
library(dplyr)
library(tidyr)
library(broom)
# For pretty plots
library(ggplot2)
theme_set(hrbrthemes::theme_ipsum())
# For nice tables
library(knitr)
# For intra-class correlations
library(irr)
# For Bayesian modeling
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# For some useful utility functions
library(hBayesDM)

library(purrr)


## ------------------------------------------------------------------------
# Data path and individual subject file names
data_path <- here::here("Data", "RawData", "Study1-Stroop")
files_t1 <- list.files(data_path, pattern = "*1.csv")
files_t2 <- list.files(data_path, pattern = "*2.csv")

# Create long-format stroop task data including all subjects
long_stroop <- map_dfr(1:47, function(i) {
  # For time 1
  tmp_t1 <- read.csv(file.path(data_path, files_t1[i]), header = F) %>%
    mutate(subj_num = i,
           time = 1)
  # For time 2 (about 3 weeks apart)
  tmp_t2 <- read.csv(file.path(data_path, files_t2[i]), header = F) %>%
    mutate(subj_num = i,
           time = 2)
  # Condition (0 = congruent, 1=neutral, 2=incongruent), 
  # Correct (1) or incorrect (0), 
  # Reaction time is in seconds
  names(tmp_t1)[1:6] <- names(tmp_t2)[1:6] <- c("Block", "Trial", "Unused", 
                                                "Condition", "Correct", "RT")
  rbind(tmp_t1, tmp_t2)
})

# Compute Stroop effect for each person at each time (equation 1)
sum_stroop <- long_stroop %>%
  group_by(subj_num, time) %>%
  summarize(stroop_eff = mean(RT[Condition==2]) - mean(RT[Condition==0]))

# Peak at the data
kable(head(sum_stroop), digits = 3)

## ------------------------------------------------------------------------
# Means and SDs of Stroop effect at each timepoint
sum_stroop %>%
  ungroup() %>%
  group_by(time) %>%
  summarize(N = n(),
            Mean = round(mean(stroop_eff), 3),
            SD = round(sd(stroop_eff), 3)) %>%
  kable(digits = 3)

## ------------------------------------------------------------------------
sum_stroop %>% 
  group_by(time) %>% 
  group_modify(~t.test(.$stroop_eff) %>% tidy) 

## ------------------------------------------------------------------------
sum_stroop %>% 
  {t.test(.$stroop_eff) %>% tidy}


## ------------------------------------------------------------------------
# Format for test-retest analysis
stroop_unpooled <- sum_stroop %>%
  ungroup() %>%
  mutate(time = ifelse(time==1, "Stroop_T1", "Stroop_T2"),
         pooled = "No", 
         Replication = "Sample Mean") %>% # these "pooled" and Replication variables will come in later
  spread(key = time, value = stroop_eff)

# Intraclass correlation of stroop effect at time 1 and 2
stroop_unpooled %>%
  select(Stroop_T1, Stroop_T2) %>%
  {irr::icc(., model = "twoway", type = "agreement", unit = "average")[c(1, 7:15)]} %>%
  as.data.frame() 

## ------------------------------------------------------------------------
stroop_unpooled %>%
  select(Stroop_T1, Stroop_T2) %>%
  irr::icc(., model = "twoway", type = "agreement", unit = "average")

## ------------------------------------------------------------------------
stroop_unpooled %>%
  select(Stroop_T1, Stroop_T2) %>%
  {cor.test(.$Stroop_T1, .$Stroop_T2)}


## ------------------------------------------------------------------------
stroop_unpooled %>% 
  ggplot(aes(Stroop_T1, Stroop_T2))+
  geom_point()+
  geom_abline(lty = 2)+
  coord_equal()+
  NULL


## ------------------------------------------------------------------------
n_subj <- long_stroop$subj_num %>% unique %>% length
n_cond <- 2
n_time <- 2
T_max <- 240

RT <- array(NA, dim = c(n_subj, n_cond, n_time, T_max))
for (i in 1:n_subj){
  RT[i, 1, 1,] = with(long_stroop, RT[subj_num==i & Condition==0 & time==1])
  # RTs for incongruent condition at time 1
  RT[i, 2, 1,] = with(long_stroop, RT[subj_num==i & Condition==2 & time==1])
  # RTs for congruent condition at time 2
  RT[i, 1, 2,] = with(long_stroop, RT[subj_num==i & Condition==0 & time==2])
  # RTs for incongruent condition at time 2
  RT[i, 2, 2,] = with(long_stroop, RT[subj_num==i & Condition==2 & time==2])
}

## ------------------------------------------------------------------------
# Stan-ready data list
stan_dat <- list(N      = n_subj,
                 N_cond = n_cond,
                 N_time = n_time,
                 T_max  = T_max,
                 RT     = RT)

## ------------------------------------------------------------------------
stroop_m1 <- stan_model("model1.stan")


## ------------------------------------------------------------------------
fit_m1 <- sampling(
  stroop_m1,
  data = stan_dat,
  iter = 2000,
  warmup = 500,
  chains = 3,
  cores = 3,
  seed = 2
)



## ------------------------------------------------------------------------
pairs(fit_m1, pars = c("mu_beta_con", "sigma_con"))

## ------------------------------------------------------------------------
stroop_m2 <- stan_model("model2.stan")

## ------------------------------------------------------------------------
fit_m2 <- sampling(stroop_m2,
                   data = stan_dat,
                   iter = 2000,
                   warmup = 500,
                   chains = 3,
                   cores = 3,
                   seed = 2
                   )


## ------------------------------------------------------------------------
pairs(fit_m2, pars = c("mu_beta_con", "sigma_con"))

## ------------------------------------------------------------------------
traceplot(fit_m1, c("mu_beta_con", "mu_beta_delta","R[1,2]"))

## ------------------------------------------------------------------------
traceplot(fit_m2, c("mu_beta_con", "mu_beta_delta","R[1,2]"))

## ------------------------------------------------------------------------
stan_rhat(fit_m1, bins = 30)
stan_rhat(fit_m2, bins = 30)


## ------------------------------------------------------------------------
pars <- extract(fit_m1)
qplot(pars$R[,1,2], geom = "density")


## ------------------------------------------------------------------------
pars <- extract(fit_m2)
qplot(pars$R[,1,2], geom = "density")


## ------------------------------------------------------------------------
pars <- extract(fit_m1)
estimate_mode(pars$R[,1,2])


## ------------------------------------------------------------------------
# Extracting posterior modes of individual-level Stroop effect estimates
stroop_pooled <- apply(pars$beta_delta, c(2,3), mean) %>%
  as.data.frame() %>%
  rename(Stroop_T1 = V1,
         Stroop_T2 = V2) %>%
  mutate(subj_num = row_number(),
         pooled = "Yes")
# Pooled correlation
pooled_cor <- with(stroop_pooled, round(cor(Stroop_T1[pooled=="Yes"], Stroop_T2[pooled=="Yes"]), 2))

# My favorite visualization of all time
bind_rows(stroop_unpooled, stroop_pooled) %>%
  mutate(subj_num = as.factor(subj_num)) %>%
  ggplot(aes(x = Stroop_T1, y = Stroop_T2)) +
  ggtitle(paste0("Posterior Means r = ", pooled_cor)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "black", size = 1) +
  stat_ellipse(geom="polygon", type="norm", level=1/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=2/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=3/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=4/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=5/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=6/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=7/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=8/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=9/10, size=0, alpha=1/10, fill="gray") +
  stat_ellipse(geom="polygon", type="norm", level=.99, size=0, alpha=1/10, fill="gray") +
  geom_line(aes(group = subj_num), size = 1/4) +
  geom_point(aes(group = subj_num, color = pooled)) +
  scale_color_manual("Pooled?",
                     values = c("#990000", "#fee8c8")) +
  theme_minimal(base_size = 20) +
  theme(panel.grid = element_blank())

