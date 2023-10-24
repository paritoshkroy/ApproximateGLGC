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
load(paste0(fpath,"./TemperatureDataAnalysis/SelectedStatelliteTemps.rda"))
head(selected.sat.temps)
table(is.na(selected.sat.temps$MaskTemp))  # FALSE are the locations to be used for modeling
nsite <- nrow(selected.sat.temps); nsite

####################################################################################
## Changing the spatial domain
####################################################################################
apply(selected.sat.temps[,c("Lon","Lat")], 2, range)  # range of the spatial domain
selected.sat.temps <- selected.sat.temps %>% 
  mutate(scaledLon = Lon - mean(range(Lon))) %>%
  mutate(scaledLat = Lat - mean(range(Lat)))
apply(selected.sat.temps[,c("scaledLon","scaledLat")], 2, range)

#####################################################################
# Preparing model objects
#####################################################################
X <- selected.sat.temps[,c("scaledLon","scaledLat")]; str(X)
y <- selected.sat.temps$TrueTemp
nsize <- selected.sat.temps %>% filter(!is.na(MaskTemp)) %>% nrow(); nsize
psize <- selected.sat.temps %>% filter(is.na(MaskTemp)) %>% nrow(); psize
obsY <- selected.sat.temps %>% filter(!is.na(MaskTemp)) %>% .$TrueTemp; str(obsY)
prdY <- selected.sat.temps %>% filter(is.na(MaskTemp)) %>% .$TrueTemp; str(prdY)
obsCoords <- selected.sat.temps %>% filter(!is.na(MaskTemp)) %>% select(scaledLon,scaledLat) %>% as.matrix() %>% unname(); str(obsCoords)
prdCoords <- selected.sat.temps %>% filter(is.na(MaskTemp)) %>% select(scaledLon,scaledLat) %>% as.matrix() %>% unname(); str(prdCoords)
obsX <- cbind(1,obsCoords); str(obsX)
prdX <- cbind(1,prdCoords); str(prdX)
################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
m1 <- 22; m2 <- 22; mstar <- m1*m2
xyRanges <- apply(selected.sat.temps[,c("scaledLon","scaledLat")], 2, range); xyRanges
Lstar <- apply(xyRanges, 2, max); Lstar
c <- c(1.5,1.5)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)
#############################################################################
# Prior elicitation
#############################################################################
obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit
rm(obsDistVec)
rm(obsDistMat)

library(nleqslv)
ab <- nleqslv(c(3,1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), from = 0, to = uLimit)

P <- 2
mu_beta <- c(mean(obsY),rep(0,P))
V_beta <- diag(c(2.5*var(obsY),rep(1,P)))
input <- list(N = nsize, M = mstar, P = 2, y = obsY, X = obsX, coords = obsCoords, L = L, lambda = lambda, mu_beta = mu_beta, V_beta = V_beta, a = ab[1], b = ab[2])
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/HSGP.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          adapt_delta = 0.99,
                          max_treedepth = 10,
                          step_size = 0.25,
                          init = 1)
elapsed_time <- cmdstan_fit$time()
elapsed_time
elapsed_time$total/3600

cmdstan_fit$cmdstan_diagnose()
sampler_diag <- cmdstan_fit$sampler_diagnostics(format = "df")
str(sampler_diag)

## Posterior summaries
pars <- c("beta[1]","beta[2]","beta[3]","sigma","ell","tau")
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
post_omega <- as_tibble(draws_df) %>% select(starts_with("omega[")) %>% as.matrix() %>% unname(); str(post_omega)
eigenfunction_compute <- function(x, L, lambda) { 
  apply(sqrt(1/L) * sin(sqrt(lambda) %*% diag(x + L)), 1, prod)
}
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
post_z <- t(sapply(1:size_post_samples, function(l) obsH %*% post_omega[l,])); str(post_z)

z_summary <- tibble(post.mean = apply(post_z, 2, mean),
                    post.sd = apply(post_z, 2, sd),
                    post.q2.5 = apply(post_z, 2, quantile2.5),
                    post.q50 = apply(post_z, 2, quantile50),
                    post.q97.5 = apply(post_z, 2, quantile97.5))
z_summary

save(elapsed_time, fixed_summary, draws_df, z_summary, file = paste0(fpath,"TemperatureDataAnalysis/HSGP_Temps.RData"))

##################################################################
## Independent prediction at each predictions sites
##################################################################

### Random effect z1 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_zpred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega[l,])); str(post_zpred)

zpred_summary <- tibble(
  post.mean = apply(post_zpred, 2, mean),
  post.sd = apply(post_zpred, 2, sd),
  post.q2.5 = apply(post_zpred, 2, quantile2.5),
  post.q50 = apply(post_zpred, 2, quantile50),
  post.q97.5 = apply(post_zpred, 2, quantile97.5))
head(zpred_summary)

### Random effect z at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_zpred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega[l,])); str(post_zpred)

zpred_summary <- tibble(
  post.mean = apply(post_zpred, 2, mean),
  post.sd = apply(post_zpred, 2, sd),
  post.q2.5 = apply(post_zpred, 2, quantile2.5),
  post.q50 = apply(post_zpred, 2, quantile50),
  post.q97.5 = apply(post_zpred, 2, quantile97.5))
head(zpred_summary)

## Compute the means
post_beta <- as_tibble(draws_df) %>% select(starts_with("beta[")) %>% as.matrix() %>% unname(); str(post_beta)
str(post_z)

str(obsX)
str(post_beta)
l <- 1
obsXbeta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_beta[l,])); str(obsXbeta)
prdXbeta <- t(sapply(1:size_post_samples, function(l) prdX %*% post_beta[l,])); str(prdXbeta)

post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
str(post_zpred)
ypred_draws <- t(sapply(1:size_post_samples, function(l) prdXbeta[l,] + post_zpred[l,] + rnorm(n = psize, mean = 0, sd = post_tau[l])))
str(ypred_draws)

pred_summary <- tibble(
  post.mean = apply(ypred_draws, 2, mean),
  post.sd = apply(ypred_draws, 2, sd),
  post.q2.5 = apply(ypred_draws, 2, quantile2.5),
  post.q50 = apply(ypred_draws, 2, quantile50),
  post.q97.5 = apply(ypred_draws, 2, quantile97.5),
  y = prdY)
head(pred_summary)

## Computation for scoring rules
## In the object PrdY there are 160 observations that are missing, to compute the scoring rules we ignore these cases
id_full_missing <- which(is.na(prdY)); str(id_full_missing)
library(scoringRules)
ES <- es_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,]); ES
VS0.25 <- vs_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,], p = 0.25); VS0.25
logs <- mean(logs_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,])); logs
CRPS <- mean(crps_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,])); CRPS

scores_df <- pred_summary %>% filter(!is.na(y)) %>%
  mutate(intervals = scoringutils::interval_score(true_values = y, lower = post.q2.5, upper = post.q50, interval_range = 0.95)) %>%
  mutate(btw = between(y,post.q2.5, post.q97.5)) %>%
  mutate(error = y - post.q50) %>%
  summarise(MAE = sqrt(mean(abs(error))), RMSE = sqrt(mean(error^2)), CVG = mean(btw),
            IS = mean(intervals)) %>%
  mutate(ES = ES, VS0.25 = VS0.25, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "HSGP") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,VS0.25,logs,`Elapsed Time`)
scores_df

save(elapsed_time, fixed_summary, draws_df, z_summary, pred_summary, scores_df, file = paste0(fpath,"TemperatureDataAnalysis/HSGP_Temps.RData"))

