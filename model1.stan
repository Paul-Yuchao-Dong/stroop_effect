data {
  int<lower=0> N; // # of subjects
  int<lower=0> N_cond; // # of conditions
  int<lower=0> N_time; // # of times of remeasurement
  int<lower=0> T_max; // # of trials on the same subject
  real RT[N, N_cond, N_time, T_max]; // reaction time
}

parameters {
  // correlation matrix for stroop effect
  corr_matrix[N_time] R;

  // SDs for group-level parameters
  vector<lower=0>[N_time] sigma_con;
  vector<lower=0>[N_time] sigma_delta;
  
  // means for group-level parameters
  vector[N_time] mu_beta_con;
  vector[N_time] mu_beta_delta;

  // SDs for normal model on RTs
  vector<lower=0>[N_cond] sigma_RT;

  // individual level parameters
  vector[N_time] beta_con[N];
  vector[N_time] beta_delta[N];
}

transformed parameters {
  // Construct covariance matrix from SDs and correlation matrix
  cov_matrix[N_time] S;
  S = quad_form_diag(R, sigma_delta);
}

model {
  // Priors on group-level SDs and correlation matrix
  R ~ lkj_corr(1);
  sigma_delta ~ normal(0, 1);
  sigma_con ~ normal(0, 1);
  sigma_RT ~ normal(0, 1);

  // Priors on group-level means
  mu_beta_con ~ normal(0,1);
  mu_beta_delta ~ normal(0,1);
  
  // Priors on individual parameters
  for (i in 1:N_time){
    beta_con[i] ~ normal(mu_beta_con[i], sigma_con[i]);
  }
  beta_delta ~ multi_normal(mu_beta_delta, S);
  
  // For each subject
  for (i in 1:N){
    // Congruent at time 1
    RT[i,1,1,:] ~ normal(beta_con[i,1], sigma_RT[1]);
    // incongruent at time 1
    RT[i,1,2,:] ~ normal(beta_con[i,1] + beta_delta[i,1], sigma_RT[2]);
    // Congruent at time 2
    RT[i,2,1,:] ~ normal(beta_con[i,2], sigma_RT[1]);
    // Incongruent at time 2
    RT[i,2,2,:] ~ normal(beta_con[i,2] + beta_delta[i,2], sigma_RT[2]);
  }
}

