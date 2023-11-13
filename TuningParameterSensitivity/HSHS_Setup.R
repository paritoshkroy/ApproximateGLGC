rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)
###########################################################################
# Local PC
###########################################################################
node <- 22
fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"
##########################################################################
# ARC Preparation
##########################################################################
fpath <- "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/" #@ARC
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
node <- as.numeric(args[1])-100 ### specify correct node here
cat("The seed used to be ", node, "\n")
##########################################################################
# Setup for the simulation study
##########################################################################
vector_lscale <- seq(0.10,0.50,l=5); vector_lscale
vector_c <- round(pmax(4.5*vector_lscale,1.5),2); vector_c

vector_m1 <- round(3.42*vector_c/vector_lscale,0); vector_m1
vector_m1 <- rep(22,length(vector_c)); vector_m1
setup1 <- tibble(lscale = vector_lscale, c = vector_c, m = vector_m1)
setup1

vector_m2 <- round(3.42*vector_c/vector_lscale,0); vector_m2
vector_m2 <- rep(32,length(vector_c)); vector_m2
setup2 <- tibble(lscale = vector_lscale, c = vector_c, m = vector_m2)
setup2

vector_m3 <- round(3.42*vector_c/vector_lscale,0); vector_m3
vector_m3 <- rep(44,length(vector_c)); vector_m3
setup3 <- tibble(lscale = vector_lscale, c = vector_c, m = vector_m3)
setup3

vector_m4 <- round(3.42*vector_c/vector_lscale,0); vector_m4
vector_m4 <- rep(51,length(vector_c)); vector_m4
setup4 <- tibble(lscale = vector_lscale, c = vector_c, m = vector_m4)
setup4

vector_m5 <- round(3.42*vector_c/vector_lscale,0); vector_m5
vector_m5 <- rep(58,length(vector_c)); vector_m4
setup5 <- tibble(lscale = vector_lscale, c = vector_c, m = vector_m5)
setup5

setup <- rbind(setup1,setup2,setup3,setup4,setup5) %>% distinct()
setup %>% print(n = nrow(setup))

##########################################################################
# Data generation
##########################################################################

this_setup <- setup[node,]; this_setup
lscale1 <- as.numeric(this_setup[,"lscale"]); lscale1
lscale2 <- as.numeric(this_setup[,"lscale"]); lscale2
c <- as.numeric(this_setup[,"c"]); c
m1 <- as.numeric(this_setup[,"m"]); m1
m2 <- as.numeric(this_setup[,"m"]); m2

source(paste0(fpath,"Rutilities/utility_functions.R"))
source(paste0(fpath,"TuningParameterSensitivity/data_generation.R"))

#######################################################################
## partition as observed and predicted
#######################################################################
obsCoords <- coords[idSampled,]
prdCoords <- coords[-idSampled,]
obsY <- y[idSampled]
prdY <- y[-idSampled]
obsX <- X[idSampled,]
prdX <- X[-idSampled,]
obsZ1 <- z1[idSampled]
obsZ2 <- z2[idSampled]
prdZ1 <- z1[-idSampled]
prdZ2 <- z2[-idSampled]

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec);obsMaxDist
obsMedDist <- median(obsDistVec);obsMedDist
obsMinDist <- min(obsDistVec);obsMinDist
rm(obsDistMat)
################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)
mstar <- m1*m2; mstar
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(c,c); c
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)

## Prior elicitation
lLimit <- as.numeric(quantile(obsDistVec, prob = 0.01)); lLimit
uLimit <- as.numeric(quantile(obsDistVec, prob = 0.50)); uLimit

## Inverse Gamma for length scale
library(nleqslv)
ab <- nleqslv(c(3,0.1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
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
input <- list(N = nsize, M = mstar, P = P, y = obsY, X = obsX, coords = obsCoords, L = L, lambda = lambda, mu_theta = mu_theta, V_theta = V_theta, a = ab[1], b = ab[2], lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, positive_skewness = 1)
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
post_omega1 <- as_tibble(draws_df) %>% select(starts_with("omega1[")) %>% as.matrix() %>% unname(); str(post_omega1)
eigenfunction_compute <- function(x, L, lambda) { 
  apply(sqrt(1/L) * sin(sqrt(lambda) %*% diag(x + L)), 1, prod)
}
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
post_z1 <- t(sapply(1:size_post_samples, function(l) obsH %*% post_omega1[l,])); str(post_z1)

z1_summary <- tibble(z1 = obsZ1,
                     post.mean = apply(post_z1, 2, mean),
                     post.sd = apply(post_z1, 2, sd),
                     post.q2.5 = apply(post_z1, 2, quantile2.5),
                     post.q50 = apply(post_z1, 2, quantile50),
                     post.q97.5 = apply(post_z1, 2, quantile97.5))
z1_summary
z1_summary %>% mutate(btw = between(z1, post.q2.5,post.q97.5)) %>% .$btw %>% mean()

## Recovery of random effect z2
size_post_samples <- nrow(draws_df); size_post_samples
post_omega2 <- as_tibble(draws_df) %>% select(starts_with("omega2[")) %>% as.matrix() %>% unname(); str(post_omega2)
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
post_z2 <- t(sapply(1:size_post_samples, function(l) obsH %*% post_omega2[l,])); str(post_z2)

z2_summary <- tibble(z2 = obsZ2,
                     post.mean = apply(post_z2, 2, mean),
                     post.sd = apply(post_z2, 2, sd),
                     post.q2.5 = apply(post_z2, 2, quantile2.5),
                     post.q50 = apply(post_z2, 2, quantile50),
                     post.q97.5 = apply(post_z2, 2, quantile97.5))
z2_summary
z2_summary %>% mutate(btw = between(z2, post.q2.5,post.q97.5)) %>% .$btw %>% mean()

save(elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, file = paste0(fpath,"TuningParameterSensitivity/HSHS_Setup",node,".RData"))

##################################################################
## Independent prediction at each predictions sites
##################################################################

### Random effect z1 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z1pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega1[l,])); str(post_z1pred)

z1pred_summary <- tibble(
  z1 = z1[-idSampled],
  post.mean = apply(post_z1pred, 2, mean),
  post.sd = apply(post_z1pred, 2, sd),
  post.q2.5 = apply(post_z1pred, 2, quantile2.5),
  post.q50 = apply(post_z1pred, 2, quantile50),
  post.q97.5 = apply(post_z1pred, 2, quantile97.5))
head(z1pred_summary)
mean(z1pred_summary[,"z1"] > z1pred_summary[,"post.q2.5"] & z1pred_summary[,"z1"] < z1pred_summary[,"post.q97.5"])

### Random effect z2 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z2pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega2[l,])); str(post_z2pred)

z2pred_summary <- tibble(
  z2 = z2[-idSampled],
  post.mean = apply(post_z2pred, 2, mean),
  post.sd = apply(post_z2pred, 2, sd),
  post.q2.5 = apply(post_z2pred, 2, quantile2.5),
  post.q50 = apply(post_z2pred, 2, quantile50),
  post.q97.5 = apply(post_z2pred, 2, quantile97.5))
head(z2pred_summary)
mean(z2pred_summary[,"z2"] > z2pred_summary[,"post.q2.5"] & z2pred_summary[,"z2"] < z2pred_summary[,"post.q97.5"])

## Compute the means
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
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
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = paste0("HSHS_Setup",node)) %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(node, c, m1, m2, elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, pred_summary, scores_df, file = paste0(fpath,"TuningParameterSensitivity/HSHS_Setup",node,".RData"))

