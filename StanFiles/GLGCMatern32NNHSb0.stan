functions {
  
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

  // 
  
  vector spdMatern32(vector lambda1, vector lambda2, real sigmasq, real lscale, int M) {
  vector[M] matern32spd;
  for(i in 1:M){
    matern32spd[i] = sigmasq * 6 * pi() * pow(sqrt(3) * inv(lscale), 3) * inv(pow(pow(sqrt(3)*inv(lscale), 2) + lambda1[i] + lambda2[i], 2.5));
  }
  return(matern32spd);
	}
	
}

data {
  int<lower=0> N;
  int<lower=0> M;
  int<lower=0> K;
  matrix[N,M] H;
  vector[M] lambda1;
  vector[M] lambda2;
  vector[N] y;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  real hypShape;
  real hypScale;
}


parameters{
  real beta;
  real<lower = 0> gamma;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> lscale1;
  real<lower = 0> lscale2;
  vector[M] noise;
}

transformed parameters{
  vector[M] omega = sqrt(spdMatern32(lambda1, lambda2, square(sigma1), lscale1, M)) .* noise;
}

model {
  beta ~ std_normal();
  gamma ~ std_normal();
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  sigma1 ~ std_normal();
  sigma2 ~ std_normal();
  tau ~ std_normal();
  noise ~ std_normal();
  y ~ vecchia_matern32(beta + gamma * exp(H * omega), square(sigma2), square(tau), lscale2, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
}



