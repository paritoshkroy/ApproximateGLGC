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
msst_df
nsite <- nrow(msst_df); nsite

#####################################################################
# Preparing model objects
#####################################################################
nsize
psize
obsY <- msst_df$temp[idSampled]; str(obsY)
prdY <- msst_df$temp[-idSampled]; str(prdY)
obsCoords <- coords[idSampled,]; str(obsCoords)
prdCoords <- coords[-idSampled,]; str(prdCoords)
obsX <- cbind(1,obsCoords); str(obsX)
prdX <- cbind(1,prdCoords); str(prdX)
str(obsCoords)
str(prdCoords)

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)
rm(obsDistMat)

################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
## For minimum m1 and m2 for the HSGP with Matern32 covariance function
quantile(obsDistVec, probs = seq(1,10,l=3)/100)
Lstar <- c(max(xRange - xMidRange), max(yRange - yMidRange)); Lstar
ell_hat <- 0.35; ell_hat
c <- pmax(1.2, 4.5*(ell_hat/Lstar)); c
L <- c*Lstar; L
m <- pmax(32,ceiling(3.42 * c/(ell_hat/Lstar))); m
m1 <- m[1]; m1
m2 <- m[2]; m2
mstar <- m1*m2; mstar

S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)
#############################################################################
# Prior elicitation
#############################################################################
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit

## Inverse Gamma for length scale
library(nleqslv)
ab <- nleqslv(c(3,2), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), 0, uLimit)

## Exponential prior for SD
lambda_sigma1 <- -log(0.01)/1; lambda_sigma1
lambda_sigma2 <- -log(0.01)/1; lambda_sigma2
lambda_tau <- -log(0.01)/1; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05

head(obsX)
P <- 3
mu_theta <- c(mean(obsY),rep(0, P-1)); mu_theta
V_theta <- diag(c(10,rep(1,P-1))); V_theta
input <- list(N = nsize, M = mstar, P = P, y = obsY, X = obsX, coords = obsCoords, L = L, lambda = lambda, mu_theta = mu_theta, V_theta = V_theta, a = ab[1], b = ab[2], lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, positive_skewness = 0)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/HSHS_GLGC_Exp.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1500,
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
post_omega1 <- as_tibble(draws_df) %>% select(starts_with("omega1[")) %>% as.matrix() %>% unname(); str(post_omega1)
eigenfunction_compute <- function(x, L, lambda) { 
  apply(sqrt(1/L) * sin(sqrt(lambda) %*% diag(x + L)), 1, prod)
}
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
post_z1 <- t(sapply(1:size_post_samples, function(l) obsH %*% post_omega1[l,])); str(post_z1)

z1_summary <- tibble(post.mean = apply(post_z1, 2, mean),
                     post.sd = apply(post_z1, 2, sd),
                     post.q2.5 = apply(post_z1, 2, quantile2.5),
                     post.q50 = apply(post_z1, 2, quantile50),
                     post.q97.5 = apply(post_z1, 2, quantile97.5))
z1_summary


## Recovery of random effect z2
size_post_samples <- nrow(draws_df); size_post_samples
post_omega2 <- as_tibble(draws_df) %>% select(starts_with("omega2[")) %>% as.matrix() %>% unname(); str(post_omega2)
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
post_z2 <- t(sapply(1:size_post_samples, function(l) obsH %*% post_omega2[l,])); str(post_z2)

z2_summary <- tibble(post.mean = apply(post_z2, 2, mean),
                     post.sd = apply(post_z2, 2, sd),
                     post.q2.5 = apply(post_z2, 2, quantile2.5),
                     post.q50 = apply(post_z2, 2, quantile50),
                     post.q97.5 = apply(post_z2, 2, quantile97.5))
z2_summary

## Recovery of z <- gamma*exp(z1) + z2
size_post_samples <- nrow(draws_df); size_post_samples
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
str(post_z1)
str(post_z2)
l <- 1
post_z <- t(sapply(1:size_post_samples, function(l) post_gamma[l]*exp(post_z1[l,]) + post_z2[l,])); str(post_z)
z_summary <- tibble(post.mean = apply(post_z, 2, mean),
                    post.sd = apply(post_z, 2, sd),
                    post.q2.5 = apply(post_z, 2, quantile2.5),
                    post.q50 = apply(post_z, 2, quantile50),
                    post.q97.5 = apply(post_z, 2, quantile97.5))
z_summary
save(elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, z_summary, post_z, file = paste0(fpath,"SSTempDataAnalysis/HSHS_GLGC_SST.RData"))

##################################################################
## Independent prediction at each predictions sites
##################################################################
## Random effect z1 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z1pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega1[l,])); str(post_z1pred)

z1pred_summary <- tibble(
  post.mean = apply(post_z1pred, 2, mean),
  post.sd = apply(post_z1pred, 2, sd),
  post.q2.5 = apply(post_z1pred, 2, quantile2.5),
  post.q50 = apply(post_z1pred, 2, quantile50),
  post.q97.5 = apply(post_z1pred, 2, quantile97.5))
z1pred_summary

### Random effect z2 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z2pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega2[l,])); str(post_z2pred)

z2pred_summary <- tibble(
  post.mean = apply(post_z2pred, 2, mean),
  post.sd = apply(post_z2pred, 2, sd),
  post.q2.5 = apply(post_z2pred, 2, quantile2.5),
  post.q50 = apply(post_z2pred, 2, quantile50),
  post.q97.5 = apply(post_z2pred, 2, quantile97.5))
z2pred_summary

## Compute the means
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(post_z)
str(post_z1)
str(post_z2)

str(obsX)
str(post_theta)

obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
prdXtheta <- t(sapply(1:size_post_samples, function(l) prdX %*% post_theta[l,])); str(prdXtheta)

post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
str(post_z1pred)
str(post_z2pred)
ypred_draws <- t(sapply(1:size_post_samples, function(l) prdXtheta[l,] + post_gamma[l] * exp(post_z1pred[l,]) + post_z2pred[l,] + rnorm(n = psize, mean = 0, sd = post_tau[l])))
str(ypred_draws)

pred_summary <- tibble(
  y = prdY,
  post.mean = apply(ypred_draws, 2, mean),
  post.sd = apply(ypred_draws, 2, sd),
  post.q2.5 = apply(ypred_draws, 2, quantile2.5),
  post.q50 = apply(ypred_draws, 2, quantile50),
  post.q97.5 = apply(ypred_draws, 2, quantile97.5))
pred_summary

## Computation for scoring rules
library(scoringRules)
ES <- es_sample(y = prdY, dat = t(ypred_draws)); ES
logs <- mean(logs_sample(y = prdY, dat = t(ypred_draws))); logs
CRPS <- mean(crps_sample(y = prdY, dat = t(ypred_draws))); CRPS

scores_df <- pred_summary %>% 
  mutate(intervals = scoringutils::interval_score(true_values = y, lower = post.q2.5, upper = post.q50, interval_range = 0.95)) %>%
  mutate(btw = between(y,post.q2.5, post.q97.5)) %>%
  mutate(error = y - post.q50) %>%
  summarise(MAE = sqrt(mean(abs(error))), RMSE = sqrt(mean(error^2)), CVG = mean(btw),
            IS = mean(intervals)) %>%
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "HSHS_GLGC") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(m1, m2, mstar, elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, z_summary, post_z, pred_summary, scores_df, file = paste0(fpath,"SSTempDataAnalysis/HSHS_GLGC_SST.RData"))

