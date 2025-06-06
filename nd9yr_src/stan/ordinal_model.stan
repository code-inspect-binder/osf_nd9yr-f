// Ordinal regression model
// 
// pushforward prior on latent cut points thanks to Michael Betancourt
// https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
//
// Matteo Lisi
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
  int<lower=1> K;             // Number of ordinal categories
  int<lower=1, upper=K> y[N]; // Observed responses
  
  vector[N] X;                // Continuous predictor (log-M-ratio)
  
  real<lower=-0.5,upper=0.5> gender[N]; // sum contrast for gender
  real<lower=-0.5,upper=0.5> covid_affected[N];
  
  vector[N] age;
  vector[N] d_prime;
  
  int<lower=0> n_grade;
  int<lower=0> n_region;
  int<lower=0> n_edu;
  int<lower=0> n_vote2019;
  int<lower=0> n_voteEUref;
  int<lower=0> n_marital_stat;
  int<lower=0> n_income;
  
  int<lower=0,upper=n_grade> grade[N];
  int<lower=0,upper=n_region> region[N];
  int<lower=0,upper=n_edu> edu[N];
  int<lower=0,upper=n_vote2019> vote2019[N];
  int<lower=0,upper=n_voteEUref> voteEUref[N];
  int<lower=0,upper=n_marital_stat> marital_stat[N];
  int<lower=0,upper=n_income> income[N];
}

transformed data{
 vector[N] age_c;
 vector[N] X_c;
 vector[N] d_prime_c;
 
 // center continuous predictors
 // and divide by 2 SD to put on similar scale as dummy variables
 age_c = (age - mean(age))/(2*sd(age));
 X_c = (X - mean(X))/(2*sd(X));
 d_prime_c = (d_prime - mean(d_prime))/(2*sd(d_prime));
}

parameters {
  real b_d;  
  real b_X;  
  real b_gender;
  real b_age;
  real b_age2;
  real b_grade[n_grade];
  real b_region[n_region];
  real b_edu[n_edu];
  real b_vote2019[n_vote2019];
  real b_voteEUref[n_voteEUref];
  real b_marital_stat[n_marital_stat];
  real b_income[n_income];
  real b_covid_affected;
  
  real<lower=0> sigma_grade;
  real<lower=0> sigma_region;
  real<lower=0> sigma_edu;
  real<lower=0> sigma_vote2019;
  real<lower=0> sigma_voteEUref;
  real<lower=0> sigma_marital_stat;
  real<lower=0> sigma_income;
  
  ordered[K - 1] c;     // cut points
}

transformed parameters {    
  vector[N] gamma;     // latent 'affinity'
  
  for (i in 1:N){
    gamma[i] = b_d * d_prime_c[i] +
               b_X * X_c[i] +
               b_gender * gender[i] +
               b_covid_affected*covid_affected[i] +
               b_age * age_c[i] + b_age2 * age_c[i]^2 +
               b_grade[grade[i]]*sigma_grade + 
               b_region[region[i]]*sigma_region + 
               b_edu[edu[i]]*sigma_edu + 
               b_vote2019[vote2019[i]]*sigma_vote2019 + 
               b_voteEUref[voteEUref[i]]*sigma_voteEUref+ 
               b_marital_stat[marital_stat[i]]*sigma_marital_stat+ 
               b_income[income[i]]*sigma_income;
  }
}

model {
  // priors
  b_X ~ normal(0, 1);
  b_d ~ normal(0, 1);
  
  c ~ induced_dirichlet(rep_vector(1, K), 0);
  
  b_age ~ normal(0, 1);
  b_age2 ~ normal(0, 1);
  b_covid_affected ~ normal(0,1);  
  b_grade ~ normal(0, 1);
  b_region ~ normal(0, 1);
  b_edu ~ normal(0, 1);
  b_vote2019 ~ normal(0, 1);
  b_voteEUref ~ normal(0, 1);
  b_marital_stat ~ normal(0, 1);
  b_income ~ normal(0, 1);
  
  // Half-cauchy priors for SDs
  sigma_grade ~ cauchy(0,1);
  sigma_region ~ cauchy(0,1);
  sigma_edu ~ cauchy(0,1);
  sigma_vote2019 ~ cauchy(0,1);
  sigma_voteEUref ~ cauchy(0,1);
  sigma_marital_stat ~ cauchy(0,1);
  sigma_income ~ cauchy(0,1);

  // likelihood
  y ~ ordered_logistic(gamma, c);
}

