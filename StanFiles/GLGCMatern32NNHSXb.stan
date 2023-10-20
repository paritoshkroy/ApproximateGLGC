functions {
  
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
  matrix[N, P + 1] X;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  matrix[N,2] coords;
  vector[2] L;
  matrix[M,2] lambda;
  vector[P + 1] mu_beta;
  matrix[P + 1, P + 1] V_beta;
  real hypShape;
  real hypScale;
}


transformed data {
  cholesky_factor_cov[P + 1] chol_V_beta;
  matrix[N,M] H;
  for(i in 1:M){
    H[,i] = eigenfunction(L, to_vector(lambda[i,]), coords);
  }
  chol_V_beta = cholesky_decompose(V_beta);
}

  
parameters{
  vector[P + 1] beta_std;
  real<lower = 0> gamma;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> lscale1;
  real<lower = 0> lscale2;
  vector[M] noise;
}

transformed parameters{
  vector[P + 1] beta;
  // implies : beta ~ multi_normal_cholesky(mu_beta, chol_V_beta);
  beta = mu_beta + chol_V_beta * beta_std;
}

model {
  vector[M] omega = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma1), lscale1, M)) .* noise;
  beta_std ~ std_normal();
  gamma ~ std_normal();
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  sigma1 ~ std_normal();
  sigma2 ~ std_normal();
  tau ~ std_normal();
  noise ~ std_normal();
  y ~ vecchia_matern32(X * beta + gamma * exp(H * omega), square(sigma2), square(tau), lscale2, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
    vector[M] omega = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma1), lscale1, M)) .* noise;
}



