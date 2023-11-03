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
  
  array[] vector predict_nnhsglgc_rng(vector y, matrix obsX, matrix predX, array[] vector obsCoords, array[] vector pred2obsDist, array[,] int pred2obsNeiID, array[] vector beta, array[] vector z1, array[] vector z1pred, vector gamma, vector sigma2, vector lscale2, vector tau, int nsize, int psize, int postsize){
    array[postsize] vector[psize] out;
    int nprint = postsize %/% 10;
    int m = dims(pred2obsDist)[2];
    int p = dims(obsX)[2];
    for(l in 1:postsize) {
      if(l%nprint == 0) print("Starts for prediction location : ", l);
      for(i in 1:psize) {
        vector[m] res = y[pred2obsNeiID[i,1:m]] - obsX[pred2obsNeiID[i,1:m],1:p]*beta[l] - gamma[l] * exp(z1[l][pred2obsNeiID[i,1:m]]);
        vector[m] c20 = matern32_distance_vector(pred2obsDist[i], 1, lscale2[l]);
        matrix[m,m] Ch2 = cholesky_decompose(add_diag(gp_matern32_cov(obsCoords[pred2obsNeiID[i,1:m]], 1, lscale2[l]), rep_vector(square(tau[l])*inv_square(sigma2[l]), m)));
        real mu2 =  predX[i,1:p]*beta[l] + gamma[l] * exp(z1pred[l][i]) + c20'*mdivide_left_tri_upp(Ch2',mdivide_left_tri_low(Ch2, res));
        real v2 = square(sigma2[l])*(1 + square(tau[l])*inv_square(sigma2[l]) - dot_self(mdivide_left_tri_low(Ch2, c20)));
        out[l][i] = normal_rng(mu2, sqrt(v2));
      }
    }
    return out;
  }

  
    // nearest neighbor approximimation of data likelihood with matern 1/2 correlation
    real vecchia_matern12_lpdf(vector y, vector Xbeta, real sigmasq, real tausq, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int N, int K) {
    
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
              neiCorMat[j, k] = exp(-inv(lscale) * neiDistMat[(i - 1), h]);
              neiCorMat[k, j] = neiCorMat[j, k];
              }
            }
          for(j in 1:dim){
            neiCorMat[j, j] = variance_ratio_plus_1;
            }
          }
        neiCorChol = cholesky_decompose(neiCorMat);
        site2neiCor = to_vector(exp(-inv(lscale) * site2neiDist[(i - 1), 1: dim]));
        v = mdivide_left_tri_low(neiCorChol, site2neiCor);
        V[i] = variance_ratio_plus_1 - dot_self(v);
        v2 = mdivide_right_tri_low(v', neiCorChol);
        U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
        }
      V[1] = variance_ratio_plus_1;
      return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) + sum(log(V)) + N * log(sigmasq));
   }
   
   // nearest neighbor approximimation of data likelihood with matern 3/2 correlation
    real vecchia_matern32_lpdf(vector y, vector Xbeta, real sigmasq, real tausq,
                             real lscale, matrix site2neiDist, matrix neiDistMat, 
                             array[,] int neiID, int N, int K) {
                               
          vector[N] V;
          vector[N] resid = y - Xbeta;
          vector[N] U = resid;
          real variance_ratio_plus_1 = tausq*inv(sigmasq) + 1; // variance ratio plus 1
          int dim;
          int h;

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
             V[i] = variance_ratio_plus_1 - dot_self(v);
             v2 = mdivide_right_tri_low(v', neiCorChol);
             U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
          }
          V[1] = variance_ratio_plus_1;
          return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
      }


  // spectral density for Hilbert space GP with matern 1/2 correlation
  vector spdMatern12(vector lambda1, vector lambda2, real sigmasq, real lscale, int M) {
    vector[M] matern12spd;
    for(i in 1:M){
      matern12spd[i] = sigmasq * 2 * pi() * square(lscale) * inv(pow(1 + square(lscale) * (lambda1[i] + lambda2[i]), 1.5));
      }
    return(matern12spd);
    }

  // spectral density for Hilbert space GP with matern 3/2 correlation
  vector spdMatern32(vector lambda1, vector lambda2, real sigmasq, real lscale, int M) {
    vector[M] matern32spd;
    for(i in 1:M){
      matern32spd[i] = sigmasq * 6 * pi() * pow(sqrt(3) * inv(lscale), 3) * inv(pow(pow(sqrt(3)*inv(lscale), 2) + lambda1[i] + lambda2[i], 2.5));
      }
    return(matern32spd);
    }
  
  // Laplacian eigenfunction for Hilbert space GP
  vector eigenfunction(vector L, vector lambda, matrix x) {
		int nc = cols(x);
		int nr = rows(x);
		matrix[nr,nc] h;
		vector[nr] h1;
		for (i in 1:nc){
			h[,i] = inv_sqrt(L[i]) * sin(sqrt(lambda[i]) * (x[,i] + L[i]));
		}
		h1 = h[,1];
		for (i in 2:nc){
			h1 = h1 .* h[,i];
		}
		return h1;
	}

}

data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0> K;
  int<lower=0> P;
  vector[N] y;
  matrix[N,P] X;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  matrix[N,2] coords;
  vector[2] L;
  matrix[M,2] lambda;
  vector[P] mu_theta;
  matrix[P,P] V_theta;
  real<lower=0> lambda_sigma1;
  real<lower=0> lambda_sigma2;
  real<lower=0> lambda_tau;
  real<lower=0> lambda_ell1;
  real<lower=0> lambda_ell2;
  int<lower=0, upper=1> positive_skewness;
}


transformed data {
  matrix[N,M] H;
  for(i in 1:M){
    H[,i] = eigenfunction(L, to_vector(lambda[i,]), coords);
  }
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
  vector[M] noise1;
}

transformed parameters{
  real gamma = skewness * abs_gamma;
  // implies : theta ~ multi_normal_cholesky(mu_theta, chol_V_theta);
  vector[P] theta = mu_theta + chol_V_theta * theta_std;
  vector[M] omega1 = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma1), ell1, M)) .* noise1;
}

model {
  vector[N] z1 = H * omega1;
  
  theta_std ~ std_normal();
  abs_gamma ~ std_normal();
  sigma1 ~ exponential(lambda_sigma1);
  sigma2 ~ exponential(lambda_sigma2);
  tau ~ exponential(lambda_tau);
  ell1 ~ frechet(1,lambda_ell1);
  ell2 ~ frechet(1,lambda_ell2);
  noise1 ~ std_normal();
  vector[N] mu = X * theta + gamma * exp(z1);
  y ~ vecchia_matern32(mu, square(sigma2), square(tau), ell2, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
  
}



