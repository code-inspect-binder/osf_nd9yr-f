// ordinal regression model for comparison 2 conditions with subject-specific 
// random effects and fully para,eterized variance-covariance matrix.
//
// Matteo Lisi, 2022
//
// pushforward prior on latent cut points thanks to Michael Betancourt
// https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
//
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data {
  int<lower=1> N;             // Number of observations
  int<lower=1> J;             // Number of participants
  int<lower=1> K;             // Number of ordinal categories
  int<lower=1, upper=K> y[N]; // confidence ratings
  int<lower=0, upper=1> x[N]; // condition, 1->Covid-19
  int<lower=0, upper=1> correct[N]; // condition, 1->Covid-19
  int<lower=1, upper=J> id[N];
}

parameters {
  real beta[2];
  ordered[K - 1] c; 
  vector<lower=0>[2] sigma_u;        // random effects standard deviations
  cholesky_factor_corr[2] L_u;       // Choleski factor of the cov matrix
  matrix[2,J] z_u;                   // random effect matrix
}

transformed parameters {
  matrix[2,J] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
}

model {
  real gamma;
  
  beta ~ normal(0, 1);
  c ~ induced_dirichlet(rep_vector(1, K), 0);

  sigma_u ~ normal(0, 1);
  L_u ~ lkj_corr_cholesky(2); 
  to_vector(z_u) ~ normal(0,1);
  
  for (i in 1:N){
    
    gamma = (beta[1] + u[1,id[i]]) * x[i] + 
            (beta[2] + u[2,id[i]]) * correct[i];
            
    y[i] ~ ordered_logistic(gamma, c);
  }
}
