data {
  int<lower=0> N;
  int<lower=0> K;
  vector[N] y;
  array[N] vector[2] coords;
  real hypShape;
  real hypScale;
}

transformed data{
  vector[N] jitter = rep_vector(1e-9, N);
}

parameters{
  real beta;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> lscale1;
  real<lower = 0> lscale2;
  vector[N] noise;
}

transformed parameters {
  }

// https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

model {
  matrix[N,N] Sigma1;
  matrix[N,N] Sigma2;
  Sigma1 = gp_matern32_cov(coords, sigma1, lscale1);
  Sigma2 = add_diag(gp_matern32_cov(coords, sigma2, lscale2), rep_vector(square(tau), N));
  
  beta ~ normal(0,10);
  sigma1 ~ normal(0,3);
  sigma2 ~ normal(0,3);
  tau ~ normal(0,2);
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  noise ~ std_normal();
  y ~ multi_normal_cholesky(X * beta + gamma * exp(cholesky_decompose(add_diag(Sigma1, jitter)) * noise), cholesky_decompose(add_diag(Sigma2, jitter)));
}

generated quantities{
  matrix[N,N] Sigma1 = gp_matern32_cov(coords, sigma1, lscale1) + diag_matrix(rep_vector(1e-9, N));
  vector[N] z = cholesky_decompose(add_diag(Sigma1, jitter)) * noise;
}

