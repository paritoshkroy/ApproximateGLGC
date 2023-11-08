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
array[] vector predict_fullgp_rng(vector y, array[] vector obsXb, array[] vector predXb, array[] vector obsCoords, array[] vector predCoords, vector sigma, vector lscale, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for posterior sample : ", l); 
      matrix[nsize,psize] C0 = gp_matern32_cov(obsCoords, predCoords, 1, lscale[l]);
      matrix[nsize,nsize] Ch = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords, 1, lscale[l]), rep_vector(square(tau[l])*inv_square(sigma[l]),nsize)));
      vector[nsize] res = y - obsXb[l];
      for(i in 1:psize) {
        row_vector[nsize] c0 = to_row_vector(C0[1:nsize,i]);
        real m =  predXb[l][i] + c0*mdivide_left_tri_upp(Ch',mdivide_left_tri_low(Ch, res));
        real v = square(sigma[l])*(1 + square(tau[l])*inv_square(sigma[l]) - dot_self(mdivide_left_tri_low(Ch, c0')));
        out[l][i] = normal_rng(m, sqrt(v));
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
  real<lower=0> lambda_sigma2;
  real<lower=0> lambda_tau;
  real<lower=0> a;
  real<lower=0> b;
}

transformed data {
  vector[N] jitter = rep_vector(1e-9, N);
  cholesky_factor_cov[P] chol_V_theta = cholesky_decompose(V_theta);
}


parameters{
  vector[P] theta_std;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> ell2;
}

transformed parameters {
  // implies : theta ~ multi_normal_cholesky(mu_theta, chol_V_theta);
  vector[P] theta = mu_theta + chol_V_theta * theta_std;
  }


model {
  matrix[N,N] C2 = gp_matern32_cov(coords, sigma2, ell2);
  matrix[N,N] L = cholesky_decompose(add_diag(C2, rep_vector(square(tau), N) + jitter));
  
  theta_std ~ std_normal();
  sigma2 ~ exponential(lambda_sigma2);
  tau ~ exponential(lambda_tau);
  ell2 ~ inv_gamma(a,b);
  vector[N] mu = X * theta;
  y ~ multi_normal_cholesky(mu, L);
}

generated quantities{
  
}

