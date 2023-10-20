functions {
  
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
          return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
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
  vector[N] y;
  array[N - 1, K] int neiID;
  matrix[N - 1, K] site2neiDist;
  matrix[N - 1, (K * (K - 1)) %/% 2] neiDistMat;
  real gamma;
  real hypShape;
  real hypScale;
}

parameters{
  real<lower = 0.0001> sigma1;
  real<lower = 0.0001> sigma2;
  real<lower = 0.0001> lscale1;
  real<lower = 0.0001> lscale2;
  real<lower = 0.0001> tau;
  vector[N] noise;
}

model {
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  sigma1 ~ normal(0,3);
  sigma2 ~ normal(0,3);
  tau ~ normal(0,2);
  noise ~ std_normal();
  y ~ vecchia_matern32(gamma * exp(latent_nngp_matern32_stuff(noise, square(sigma1), lscale1, site2neiDist, neiDistMat, neiID, N, K)), square(sigma2), square(tau), lscale2, site2neiDist, neiDistMat, neiID, N, K);
}

generated quantities {
  vector[N] z = latent_nngp_matern32_stuff(noise, square(sigma1), lscale1, site2neiDist, neiDistMat, neiID, N, K);
}

