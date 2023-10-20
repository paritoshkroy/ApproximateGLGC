functions {
  
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
  vector[N] y;
  matrix[N,2] coords;
  vector[2] L;
  matrix[M,2] lambda;
  real mu_beta;
  real<lower=0> V_beta;
  real<lower=0> hypShape;
  real<lower=0> hypScale;
}

transformed data {
  real sqrt_V_beta;
  matrix[N,M] H;
  for(i in 1:M){
    H[,i] = eigenfunction(L, to_vector(lambda[i,]), coords);
  }
  sqrt_V_beta = sqrt(V_beta);
}


parameters{
  real beta_raw;
  real<lower = 0> gamma;
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> lscale1;
  real<lower = 0> lscale2;
  vector[M] noise1;
  vector[M] noise2;
}

transformed parameters{
  real beta;
  // implies : beta ~ normal(mu_beta, sqrt_V_beta);
  beta = mu_beta + sqrt_V_beta * beta_raw;
  
}

model {
  
  vector[M] omega1 = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma1), lscale1, M)) .* noise1;
  vector[M] omega2 = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma2), lscale2, M)) .* noise2;
  vector[N] z1 = H * omega1;
  vector[N] z2 = H * omega2;
  
  beta_raw ~ std_normal();
  gamma ~ std_normal();
  sigma1 ~ std_normal();
  sigma2 ~ std_normal();
  tau ~ std_normal();
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  noise1 ~ std_normal();
  noise2 ~ std_normal();
  y ~ normal(beta + gamma * exp(z1) + z2, tau);
}

generated quantities {
  vector[M] omega1 = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma1), lscale1, M)) .* noise1;
  vector[M] omega2 = sqrt(spdMatern32(lambda[,1], lambda[,2], square(sigma2), lscale2, M)) .* noise2;
}




