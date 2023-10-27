rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"
fpath <- "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/" #@ARC

source(paste0(fpath,"Rutilities/utility_functions.R"))
source(paste0(fpath,"ExactVsApproximateMethods/data_generation.R"))

# partition as observed and predicted
obsCoords <- coords[idSampled,]
prdCoords <- coords[-idSampled,]
obsY <- y[idSampled]
prdY <- y[-idSampled]
obsX <- X[idSampled,]
prdX <- X[-idSampled,]

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)

# Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.025); lLimit
uLimit <- quantile(obsDistVec, prob = 0.975); uLimit
rm(obsDistMat)

lambda_sigma1 <- -log(0.01)/1; lambda_sigma1
lambda_sigma2 <- -log(0.01)/1; lambda_sigma2
lambda_tau <- -log(0.01)/1; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05
lambda_ell1 <- as.numeric(-log(0.01)*lLimit); lambda_ell1
lambda_ell2 <- as.numeric(-log(0.01)*lLimit); lambda_ell2
pfrechet(q = lLimit, alpha = 1, sigma = lambda_ell2, lower.tail = TRUE) ## P(ell < lLimit) = 0.05

P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
V_theta <- diag(c(10,rep(1,P-1)))
input <- list(N = nsize, P = P, y = obsY, X = obsX, coords = obsCoords, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, lambda_ell1 = lambda_ell1, lambda_ell2 = lambda_ell2, positive_skewness = 1)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/Full_GLGC_PC.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1500,
                          iter_sampling = 500,
                          adapt_delta = 0.99,
                          max_treedepth = 12,
                          step_size = 0.25)
elapsed_time <- cmdstan_fit$time()
elapsed_time
elapsed_time$total/3600

cmdstan_fit$cmdstan_diagnose()
sampler_diag <- cmdstan_fit$sampler_diagnostics(format = "df")
str(sampler_diag)

## Posterior summaries
pars <- c(paste0("theta[",1:P,"]"),"sigma1","sigma2","ell1","ell2","tau","gamma")
pars_true_df <- tibble(variable = pars, true = c(theta,sigma1,sigma2,lscale1,lscale2,tau,gamma))
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- inner_join(pars_true_df, fit_summary)
fixed_summary %>% print(digits = 3)

## Posterior draws
draws_df <- cmdstan_fit$draws(format = "df")
draws_df

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + facet_text(size = 15)

## Recovery of random effect z1
size_post_samples <- nrow(draws_df); size_post_samples
post_ell1 <- as_tibble(draws_df) %>% .$ell1; str(post_ell1)
post_sigma1 <- as_tibble(draws_df) %>% .$sigma1; str(post_sigma1)
post_noise1 <- as_tibble(draws_df) %>% select(starts_with("noise1[")) %>% as.matrix() %>% unname(); str(post_noise1)
post_z1 <- array(0, dim = c(size_post_samples,nsize)); str(post_z1)
obsDistMat <- fields::rdist(obsCoords)
l <- 1
for(l in 1:size_post_samples){
  C1 <- matern32(d = obsDistMat, sigma = post_sigma1[l],  lscale = post_ell1[l])
  post_z1[l,] <- drop(crossprod(chol(C1),post_noise1[l,]))
}
str(post_z1)

z1_summary <- tibble(z1 = z1[idSampled],
                     post.mean = apply(post_z1, 2, mean),
                     post.sd = apply(post_z1, 2, sd),
                     post.q2.5 = apply(post_z1, 2, quantile2.5),
                     post.q50 = apply(post_z1, 2, quantile50),
                     post.q97.5 = apply(post_z1, 2, quantile97.5))
z1_summary
z1_summary %>% mutate(btw = between(z1, post.q2.5,post.q97.5)) %>% .$btw %>% mean()

save(elapsed_time, fixed_summary, draws_df, z1_summary, file = paste0(fpath,"ExactVsApproximateMethods/Full_GLGC_PC.RData"))

##################################################################
## Independent prediction at each predictions sites
##################################################################
source(paste0(fpath,"Rutilities/expose_cmdstanr_functions.R"))
exsf <- expose_cmdstanr_functions(model_path = stan_file)
args(exsf$predict_fullglgc_rng)

size_post_samples <- nrow(draws_df); size_post_samples
psize <- nrow(prdCoords); psize

post_sigma1 <- as_tibble(draws_df) %>% .$sigma1; str(post_sigma1)
post_sigma2 <- as_tibble(draws_df) %>% .$sigma2; str(post_sigma2)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell1 <- as_tibble(draws_df) %>% .$ell1; str(post_ell1)
post_ell2 <- as_tibble(draws_df) %>% .$ell2; str(post_ell2)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(post_z1)

str(obsX)
str(post_theta)
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
prdXtheta <- t(sapply(1:size_post_samples, function(l) prdX %*% post_theta[l,])); str(prdXtheta)

str(exsf$my_gp_matern32_cov(x = lapply(1:nsize, function(i) obsCoords[i,]), y = lapply(1:psize, function(i) prdCoords[i,]), sigma = 1, lscale = 1))

post_ypred <- exsf$predict_fullglgc_rng(
  y = obsY, 
  obsXb = lapply(1:size_post_samples, function(i) obsXtheta[i,]), 
  predXb = lapply(1:size_post_samples, function(i) prdXtheta[i,]), 
  obsCoords = lapply(1:nsize, function(i) obsCoords[i,]), 
  predCoords = lapply(1:psize, function(i) prdCoords[i,]), 
  z1 = lapply(1:size_post_samples, function(i) post_z1[i,]), 
  gamma = post_gamma, 
  sigma1 = post_sigma1, 
  sigma2 = post_sigma2, 
  lscale1 = post_ell1,
  lscale2 = post_ell2, 
  tau = post_tau, 
  nsize = nsize, 
  psize = psize,
  postsize = size_post_samples)

ypred_draws <- do.call(rbind,post_ypred); str(ypred_draws)
pred_summary <- tibble(
  post.mean = apply(ypred_draws, 2, mean),
  post.sd = apply(ypred_draws, 2, sd),
  post.q2.5 = apply(ypred_draws, 2, quantile2.5),
  post.q50 = apply(ypred_draws, 2, quantile50),
  post.q97.5 = apply(ypred_draws, 2, quantile97.5),
  y = prdY)
head(pred_summary)
mean(pred_summary[,"y"]>pred_summary[,"post.q2.5"] & pred_summary[,"y"]<pred_summary[,"post.q97.5"])

## Computation for scoring rules

library(scoringRules)
ES <- es_sample(y = prdY, dat = t(ypred_draws)); ES
#VS0.25 <- vs_sample(y = prdY, dat = t(ypred_draws), p = 0.25); VS0.25
logs <- mean(logs_sample(y = prdY, dat = t(ypred_draws))); logs
CRPS <- mean(crps_sample(y = prdY, dat = t(ypred_draws))); CRPS

scores_df <- pred_summary %>% 
  mutate(intervals = scoringutils::interval_score(true_values = y, lower = post.q2.5, upper = post.q50, interval_range = 0.95)) %>%
  mutate(btw = between(y,post.q2.5, post.q97.5)) %>%
  mutate(error = y - post.q50) %>%
  summarise(MAE = sqrt(mean(abs(error))), RMSE = sqrt(mean(error^2)), CVG = mean(btw),
            IS = mean(intervals)) %>%
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "Full_GLGC_PC") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(elapsed_time, fixed_summary, draws_df, z1_summary, pred_summary, scores_df, file = paste0(fpath,"ExactVsApproximateMethods/Full_GLGC_PC.RData"))

