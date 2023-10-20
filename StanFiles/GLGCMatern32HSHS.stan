functions {
  
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
  vector[N] y;
  real gamma;
  matrix[N,M] H;
  vector[M] lambda1;
  vector[M] lambda2;
  real hypShape;
  real hypScale;
}


parameters{
  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> tau;
  real<lower = 0> lscale1;
  real<lower = 0> lscale2;
  vector[M] noise1;
  vector[M] noise2;
}

transformed parameters{
  vector[M] omega1 = sqrt(spdMatern32(lambda1, lambda2, square(sigma1), lscale1, M)) .* noise1;
  vector[M] omega2 = sqrt(spdMatern32(lambda1, lambda2, square(sigma2), lscale2, M)) .* noise2;
}

model {
  sigma1 ~ normal(0,3);
  sigma2 ~ normal(0,3);
  tau ~ normal(0,2);
  lscale1 ~ inv_gamma(hypShape,hypScale);
  lscale2 ~ inv_gamma(hypShape,hypScale);
  noise1 ~ std_normal();
  noise2 ~ std_normal();
  y ~ normal(gamma * exp(H * omega1) + H * omega2, tau);
}

generated quantities {
  
}



