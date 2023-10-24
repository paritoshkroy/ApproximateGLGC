functions {
  
  matrix my_gp_matern32_cov(array[] vector x, array[] vector y, real sigma, real lscale){
  return gp_matern32_cov(x, y, sigma, lscale);
  }
  
  vector matern32_distance_vector(vector d, real sigma, real lscale){
   vector[size(d)] out;
   vector[size(d)] ds = sqrt(3)*d*inv(lscale);
   out = square(sigma) * ((1 + ds) .* exp(-ds));
   return out;
  }
  
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
  
  array[] vector predict_nnnnglgc_rng(vector y, matrix obsX, matrix predX, array[] vector obsCoords, array[] vector pred2obsDist, array[,] int pred2obsNeiID, array[] vector beta, array[] vector z1, vector gamma, vector sigma1, vector sigma2, vector lscale1, vector lscale2, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    int m = dims(pred2obsDist)[2];
    int p = dims(obsX)[2];
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for prediction location : ", l);
      for(i in 1:psize) {
        vector[m] c10 = matern32_distance_vector(pred2obsDist[i], 1, lscale1[l]);
        matrix[m,m] Ch1 = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords[pred2obsNeiID[i,1:m]], 1, lscale1[l]), rep_vector(1e-7, m)));
        real mu1 =  c10'*mdivide_left_tri_upp(Ch1',mdivide_left_tri_low(Ch1, z1[l][pred2obsNeiID[i,1:m]]));
        real v1 = square(sigma1[l])*(1 - dot_self(mdivide_left_tri_low(Ch1, c10)));
        real z10 = normal_rng(mu1, sqrt(v1));
        vector[m] res = y[pred2obsNeiID[i,1:m]] - obsX[pred2obsNeiID[i,1:m],1:p]*beta[l] - gamma[l] * exp(z1[l][pred2obsNeiID[i,1:m]]);
        vector[m] c20 = matern32_distance_vector(pred2obsDist[i], 1, lscale2[l]);
        matrix[m,m] Ch2 = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords[pred2obsNeiID[i,1:m]], 1, lscale2[l]), rep_vector(square(tau[l])*inv_square(sigma2[l]), m)));
        real mu2 =  predX[i,1:p]*beta[l] + gamma[l] * exp(z10) + c20'*mdivide_left_tri_upp(Ch2',mdivide_left_tri_low(Ch2, res));
        real v2 = square(sigma2[l])*(1 + square(tau[l])*inv_square(sigma2[l]) - dot_self(mdivide_left_tri_low(Ch2, c20)));
        out[l][i] = normal_rng(mu2, sqrt(v2));
      }
    }
    return out;
  }
  
  real vecchia_matern32_lpdf(vector y, vector Xbeta, real sigmasq, real tausq, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int N, int K) {
    
    vector[N] V;
    vector[N] resid = y - Xbeta;
    vector[N] U = resid;
    real variance_ratio_plus_1 = tausq * inv(sigmasq) + 1; // variance ratio plus 1
    int h;
    for (i in 2:N) {
      int dim = (i < (K + 1))? (i - 1) : K;
      matrix[dim, dim] neiCorMat;
      matrix[dim, dim] neiCorChol;
      vector[dim] site2neiCor;
      vector[dim] v;
      row_vector[dim] v2;
      
      if(dim == 1){
        neiCorMat[1, 1] = variance_ratio_plus_1;
        } else {
          h = 0;
          for (j in 1:(dim - 1)){
            for (k in (j + 1):dim){
              h = h + 1;
              neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
              neiCorMat[k, j] = neiCorMat[j, k];
              }
            }
            for(j in 1:dim){
              neiCorMat[j, j] = variance_ratio_plus_1;
            }
        }

        neiCorChol = cholesky_decompose(neiCorMat);
        site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)));
        v = mdivide_left_tri_low(neiCorChol, site2neiCor);
        V[i] = variance_ratio_plus_1 - dot_self(v);
        v2 = mdivide_right_tri_low(v', neiCorChol);
        U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
        }
        V[1] = variance_ratio_plus_1;
        return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) + sum(log(V)) + N * log(sigmasq));
      }
  
  vector latent_nngp_matern32_stuff(vector noise, real sigmasq, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int N, int K){
    
    
    vector[N] z;
    vector[N] z_var;
    int h;
    
    z_var[1] = sigmasq;
    z[1] = sqrt(z_var[1]) * noise[1];
  
    for (i in 2:N) {
      int dim = (i < (K + 1))? (i - 1) : K;
      matrix[dim, dim] neiCorMat;
      matrix[dim, dim] neiCorChol;
      vector[dim] site2neiCor;
      vector[dim] v;
      
      if(dim == 1){
        neiCorMat[1, 1] = 1;
      }else {
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
            neiCorMat[k, j] = neiCorMat[j, k];
          }
        }
        for(j in 1:dim){
          neiCorMat[j, j] = 1;
        }
      }
    neiCorChol = cholesky_decompose(neiCorMat);
    site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1)][1:dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1)][1:dim] * inv(lscale)));
    v = mdivide_left_tri_low(neiCorChol, site2neiCor);
    z_var[i] = sigmasq*(1-dot_self(v));
    z[i] = mdivide_right_tri_low(v', neiCorChol) * z[neiID[(i - 1), 1:dim]] + sqrt(z_var[i]) * noise[i];
  }
  return z;
}
}

data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> P;
  vector[N] y;
  matrix[N, P + 1] X;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  vector[P + 1] mu_beta;
  matrix[P + 1, P + 1] V_beta;
  real<lower=0> lambda_sigma1;
  real<lower=0> lambda_sigma2;
  real<lower=0> lambda_tau;
  real<lower=0> lambda_ell1;
  real<lower=0> lambda_ell2;
  int<lower=0, upper=1> positive_skewness;
}

transformed data {
  cholesky_factor_cov[P + 1] chol_V_beta;
  chol_V_beta = cholesky_decompose(V_beta);
  int skewness;
  if(positive_skewness==0){
    skewness = -1;
    } else {
      skewness = 1;
    }
}

parameters{
  vector[P + 1] beta_std;
  real<lower = 0> absgamma;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> ell1;
  real<lower = 0> ell2;
  real<lower = 0> tau;
  vector[N] noise;
}

transformed parameters{
  // implies : beta ~ multi_normal_cholesky(mu_beta, chol_V_beta);
  vector[P + 1] beta = mu_beta + chol_V_beta * beta_std;
  real gamma = skewness * absgamma;
}

model {
  beta_std ~ std_normal();
  absgamma ~ std_normal();
  sigma1 ~ exponential(lambda_sigma1);
  sigma2 ~ exponential(lambda_sigma2);
  tau ~ exponential(lambda_tau);
  ell1 ~ frechet(1,lambda_ell1);
  ell2 ~ frechet(1,lambda_ell2);
  noise ~ std_normal();
  y ~ vecchia_matern32(X * beta + gamma * exp(latent_nngp_matern32_stuff(noise, square(sigma1), ell1, site2neiDist, neiDistMat, neiID, N, K)), square(sigma2), square(tau), ell2, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
  
}



