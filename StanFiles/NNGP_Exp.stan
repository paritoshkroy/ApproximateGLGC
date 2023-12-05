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
  
  // fitted values using nearest neighbor approximimation of data likelihood with matern 3/2 
    vector vecchia_matern32_fitted_rng(vector y, vector mu, real sigmasq, real tausq,
                             real lscale, matrix site2neiDist, matrix neiDistMat, 
                             array[,] int neiID, int N, int K) {
                               
          vector[N] yfitted;
          vector[N] cond_mu; // conditional mean
          vector[N] V; // conditional variances = sigmasq*V
          vector[N] resid = y - mu;
          //vector[N] U = resid;
          real variance_ratio_plus_1 = tausq*inv(sigmasq) + 1; // variance ratio plus 1
          V[1] = variance_ratio_plus_1;
          int dim;
          int h;
          yfitted[1] = normal_rng(mu[1], sqrt(sigmasq*V[1]));
          
          for (i in 2:N) {
            dim = (i < (K + 1))? (i - 1) : K;
            matrix[dim, dim] neiCorMat;
              matrix[dim, dim] neiCorChol;
              vector[dim] site2neiCor;
              vector[dim] v;
              row_vector[dim] v2;

              if(dim == 1){neiCorMat[1, 1] = variance_ratio_plus_1;}
              else{
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
             V[i] = variance_ratio_plus_1 - dot_self(v); // conditional variances
             v2 = mdivide_right_tri_low(v', neiCorChol);
             cond_mu[i] = mu[i] + v2 * resid[neiID[(i - 1), 1:dim]];
             yfitted[i] = normal_rng(cond_mu[i], sqrt(sigmasq*V[i]));
          }
          return yfitted;
      }

  array[] vector predict_vecchia_rng(vector y, matrix obsX, matrix predX, array[] vector obsCoords, array[] vector pred2obsDist, array[,] int pred2obsNeiID, array[] vector beta, vector sigma, vector lscale, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    int m = dims(pred2obsDist)[2];
    int p = dims(obsX)[2];
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for prediction location : ", l);
      for(i in 1:psize) {
        vector[m] res = y[pred2obsNeiID[i,1:m]] - obsX[pred2obsNeiID[i,1:m],1:p]*beta[l];
        vector[m] c0 = matern32_distance_vector(pred2obsDist[i], 1, lscale[l]);
        matrix[m,m] Ch = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords[pred2obsNeiID[i,1:m]], 1, lscale[l]), rep_vector(square(tau[l])*inv_square(sigma[l]), m)));
        real mu =  predX[i,1:p]*beta[l] + c0'*mdivide_left_tri_upp(Ch',mdivide_left_tri_low(Ch, res));
        real v = square(sigma[l])*(1 + square(tau[l])*inv_square(sigma[l]) - dot_self(mdivide_left_tri_low(Ch, c0)));
        out[l][i] = normal_rng(mu, sqrt(v));
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
  matrix[N, P] X;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  vector[P] mu_theta;
  matrix[P, P] V_theta;
  real<lower=0> lambda_sigma;
  real<lower=0> lambda_tau;
  real a;
  real b;
}

transformed data {
  cholesky_factor_cov[P] chol_V_theta;
  chol_V_theta = cholesky_decompose(V_theta);
}

parameters{
  vector[P] theta_std;
  real<lower = 0> sigma;
  real<lower = 0> ell;
  real<lower = 0> tau;
}

transformed parameters{
  // implies : theta ~ multi_normal_cholesky(mu_theta, chol_V_theta);
  vector[P] theta = mu_theta + chol_V_theta * theta_std;
}

model {
  theta_std ~ std_normal();
  sigma ~ exponential(lambda_sigma);
  ell ~ inv_gamma(a,b);
  tau ~ exponential(lambda_tau);
  y ~ vecchia_matern32(X * theta, square(sigma), square(tau), ell, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
  
}



