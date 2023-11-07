functions{ // functions starts

matrix my_gp_matern32_cov(array[] vector x, array[] vector y, real sigma, real lscale){
  return gp_matern32_cov(x, y, sigma, lscale);
}

/* stan function for backsolve*/
vector mdivide_left_tri_upp(matrix U, vector b) { 
  int n = rows(U);
  vector[n] x = b; 
  real cs;
  array[n] int ids = sort_desc(linspaced_int_array(n, 1, n));
  for (i in ids){
    x[i] = x[i]/U[i,i];
    cs = 0;
    for (j in (i+1):n){
      cs = cs + x[j]*U[i,j];
  }
  x[i] = x[i] - cs/U[i,i];
}
return x;
}
 
/* */
array[] vector predict_fullglgc_rng(vector y, array[] vector obsXb, array[] vector predXb, array[] vector obsCoords, array[] vector predCoords, array[] vector z1obs, vector gamma, vector sigma1, vector sigma2, vector lscale1, vector lscale2, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for posterior sample : ", l); 
      matrix[nsize,psize] C10 = gp_matern32_cov(obsCoords, predCoords, 1, lscale1[l]);
      matrix[nsize,psize] C20 = gp_matern32_cov(obsCoords, predCoords, 1, lscale2[l]);
      matrix[nsize,nsize] Ch1 = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords, 1, lscale1[l]), rep_vector(1e-7,nsize)));
      matrix[nsize,nsize] Ch2 = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords, 1, lscale2[l]), rep_vector(square(tau[l])*inv_square(sigma2[l]),nsize)));
      vector[nsize] res = y - obsXb[l] - gamma[l]*exp(z1obs[l]);
      for(i in 1:psize) {
        row_vector[nsize] c10 = to_row_vector(C10[1:nsize,i]);
        real m1 =  c10*mdivide_left_tri_upp(Ch1',mdivide_left_tri_low(Ch1, z1obs[l]));
        real v1 = square(sigma1[l])*(1 - dot_self(mdivide_left_tri_low(Ch1, c10')));
        real z10 = normal_rng(m1, sqrt(v1));
        row_vector[nsize] c20 = to_row_vector(C20[1:nsize,i]);
        real m2 =  predXb[l][i] + gamma[l] * exp(z10) + c20*mdivide_left_tri_upp(Ch2',mdivide_left_tri_low(Ch2, res));
        real v2 = square(sigma2[l])*(1 + square(tau[l])*inv_square(sigma2[l]) - dot_self(mdivide_left_tri_low(Ch2, c20')));
        out[l][i] = normal_rng(m2, sqrt(v2));
      }
    }
    return out;
  }
  
} // functions ends

data {
  int<lower=0> N;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
  array[N] vector[2] coords;
  vector[P] mu_theta;
  matrix[P,P] V_theta;
  real<lower=0> lambda_sigma1;
  real<lower=0> lambda_sigma2;
  real<lower=0> lambda_tau;
  real<lower=0> a;
  real<lower=0> b;
  int<lower=0, upper=1> positive_skewness;
}

transformed data {
  vector[N] jitter = rep_vector(1e-9, N);
  cholesky_factor_cov[P] chol_V_theta = cholesky_decompose(V_theta);
  int skewness;
  if(positive_skewness==0){
    skewness = -1;
    } else {
      skewness = 1;
    }
}


parameters{
  vector[P] theta_std;
  real<lower = 0> abs_gamma;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> ell1;
  real<lower = 0> ell2;
  vector[N] noise1;
}

transformed parameters {
  real gamma = skewness * abs_gamma;
  // implies : theta ~ multi_normal_cholesky(mu_theta, chol_V_theta);
  vector[P] theta = mu_theta + chol_V_theta * theta_std;
  }


model {
  matrix[N,N] C1 = gp_matern32_cov(coords, sigma1, ell1);
  matrix[N,N] C2 = gp_matern32_cov(coords, sigma2, ell2);
  vector[N] z1 = cholesky_decompose(add_diag(C1,jitter)) * noise1;
  matrix[N,N] L = cholesky_decompose(add_diag(C2, rep_vector(square(tau), N) + jitter));
  
  theta_std ~ std_normal();
  abs_gamma ~ std_normal();
  sigma1 ~ exponential(lambda_sigma1);
  sigma2 ~ exponential(lambda_sigma2);
  tau ~ exponential(lambda_tau);
  ell1 ~ inv_gamma(a,b);
  ell2 ~ inv_gamma(a,b);
  noise1 ~ std_normal();
  vector[N] mu = X * theta + gamma * exp(z1);
  y ~ multi_normal_cholesky(mu, L);
}

generated quantities{
  
}

