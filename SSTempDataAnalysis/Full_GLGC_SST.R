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
load(paste0(fpath,"SSTempDataAnalysis/SelectedData/SSTempDataPreparation.rda"))
eastern_msst
nsite <- nrow(eastern_msst); nsite

#####################################################################
# Preparing model objects
#####################################################################
X <- cbind(1,unname(as.matrix(eastern_msst[,c("relocateLon","relocateLat")]))); str(X)
y <- eastern_msst$temp
nsize <- length(idSampled); nsize
psize <- nsite - nsize; psize
obsY <- y[idSampled]; str(obsY)
prdY <- y[-idSampled]; str(prdY)
obsCoords <- scaled.coords[idSampled,]; str(obsCoords)
prdCoords <- scaled.coords[-idSampled,]; str(prdCoords)
obsX <- X[idSampled,]; str(obsX)
prdX <- X[-idSampled,]; str(prdX)

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)
rm(obsDistMat)

## Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.05); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit

lambda_sigma1 <- -log(0.01)/1; lambda_sigma1
lambda_sigma2 <- -log(0.01)/1; lambda_sigma2
lambda_tau <- -log(0.01)/1; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05

library(nleqslv)
ab <- nleqslv(c(3,1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), 0, 1.5*uLimit)

P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
V_theta <- diag(c(10,rep(1,P-1)))

# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, P = P, y = obsY, X = obsX, coords = obsCoords, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, a = ab[1], b = ab[2], positive_skewness = 0)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/Full_GLGC_Exp.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input,
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          adapt_delta = 0.99,
                          max_treedepth = 15,
                          step_size = 0.25)
elapsed_time <- cmdstan_fit$time()
elapsed_time
elapsed_time$total/3600

cmdstan_fit$cmdstan_diagnose()
sampler_diag <- cmdstan_fit$sampler_diagnostics(format = "df")
str(sampler_diag)

## Posterior summaries
pars <- c(paste0("theta[",1:P,"]"),"sigma1","sigma2","ell1","ell2","tau","gamma")
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fit_summary %>% filter(variable %in% pars)
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

z1_summary <- tibble(post.mean = apply(post_z1, 2, mean),
                     post.sd = apply(post_z1, 2, sd),
                     post.q2.5 = apply(post_z1, 2, quantile2.5),
                     post.q50 = apply(post_z1, 2, quantile50),
                     post.q97.5 = apply(post_z1, 2, quantile97.5))
z1_summary

save(elapsed_time, fixed_summary, draws_df, z1_summary, file = paste0(fpath,"SSTempDataAnalysis/Full_GLGC_SST.RData"))

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
pred_summary
mean(pred_summary[,"y"]>pred_summary[,"post.q2.5"] & pred_summary[,"y"]<pred_summary[,"post.q97.5"])

## Computation for scoring rules
## In the object PrdY there are 160 observations that are missing, to compute the scoring rules we ignore these cases
library(scoringRules)
ES <- es_sample(y = prdY, dat = t(ypred_draws)); ES
logs <- mean(logs_sample(y = prdY, dat = t(ypred_draws))); logs
CRPS <- mean(crps_sample(y = prdY, dat = t(ypred_draws))); CRPS

scores_df <- pred_summary %>% filter(!is.na(y)) %>%
  mutate(intervals = scoringutils::interval_score(true_values = y, lower = post.q2.5, upper = post.q50, interval_range = 0.95)) %>%
  mutate(btw = between(y,post.q2.5, post.q97.5)) %>%
  mutate(error = y - post.q50) %>%
  summarise(MAE = sqrt(mean(abs(error))), RMSE = sqrt(mean(error^2)), CVG = mean(btw),
            IS = mean(intervals)) %>%
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "Full_GLGC") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(elapsed_time, fixed_summary, draws_df, z1_summary, pred_summary, scores_df, file = paste0(fpath,"SSTempDataAnalysis/Full_GLGC_SST.RData"))
