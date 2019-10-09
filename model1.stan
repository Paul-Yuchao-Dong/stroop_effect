data {
  int<lower=0> N; // # of subjects
  int<lower=0> N_cond; // # of conditions
  int<lower=0> N_time; // # of times of remeasurement
  int<lower=0> T_max; // # of trials on the same subject
  real RT[N, N_cond, N_time, T_max]; // reaction time
}

parameters {
  // Group-level correlation matrix cholesky
  cholesky_factor_corr[2] R_cholesky;

  // SDs for group-level parameters
  vector<lower=0>[2] sigma_con;
  vector<lower=0>[2] sigma_delta;
  
  // means for group-level parameters
  vector[2] mu_beta_con;
  vector[2] mu_beta_delta;

  // SDs for normal model on RTs
  vector<lower=0>[2] sigma_RT;

  // individual level parameters (before being transformed)
  matrix[N,2] beta_con_pr;
  matrix[2,N] beta_delta_pr; //order flipped here for operation?
}

transformed parameters {
  
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[2,N] beta_delta_tilde; 
  
  // Individual-level parameters
  matrix[N,2] beta_con;
  matrix[N,2] beta_delta;
  
  // Construct individual offsets
  beta_delta_tilde = diag_pre_multiply(sigma_delta, R_cholesky) * beta_delta_pr;

  for (i in 1:N){
    // Congruent at time 1
    beta_con[i, 1] = mu_beta_con[1] + sigma_con[1]*beta_con_pr[i,1];
    // Congruent at time 2
    beta_con[i, 2] = mu_beta_con[2] + sigma_con[2]*beta_con_pr[i,2];
    // Stroop effect at time 1
    beta_delta[i,1] = mu_beta_delta[1] + beta_delta_tilde[1, i];
    // Stroop effect at time 2
    beta_delta[i,2] = mu_beta_delta[2] + beta_delta_tilde[2,i];
  }
}

model {
  // Priors on cholesky factor of correlation matrix
  R_cholesky ~ lkj_corr_cholesky(1);
  
  // Priors on group-level SDs
  sigma_delta ~ cauchy(0, 1);
  sigma_con ~ cauchy(0, 1);
  sigma_RT ~ cauchy(0, 1);

  
  // Priors on individual parameters
  to_vector(beta_delta_pr) ~ normal(0,1);
  to_vector(beta_con_pr) ~ normal(0,1);
  
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

generated quantities{
  corr_matrix[2] R;
  // Reconstruct correlation matrix from cholesky factor
  R = R_cholesky * R_cholesky';
}
