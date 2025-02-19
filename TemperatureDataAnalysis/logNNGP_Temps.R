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
load(paste0(fpath,"/TemperatureDataAnalysis/SelectedData/SelectedStatelliteTemps.rda"))
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
## NNGP preparation
################################################################################
source(paste0(fpath,"Rutilities/NNMatrix.R"))
nNeighbors <- 10
neiMatInfo <- NNMatrix(coords = obsCoords, n.neighbors = nNeighbors, n.omp.threads = 2)
str(neiMatInfo)
obsY <- obsY[neiMatInfo$ord] # ordered the data following neighborhood settings
obsX <- obsX[neiMatInfo$ord,] # ordered the data following neighborhood settings
obsCoords <- obsCoords[neiMatInfo$ord,] # ordered the data following neighborhood settings

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)
lLimit <- quantile(obsDistVec, prob = 0.01)/2.75; lLimit
uLimit <- quantile(obsDistVec, prob = 0.50)/2.75; uLimit
rm(obsDistMat)

#############################################################################
# Prior elicitation
#############################################################################
library(nleqslv)
ab <- nleqslv(c(5,0.1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), from = 0, to = uLimit)

lambda_sigma <- -log(0.01)/0.01; lambda_sigma
lambda_tau <- -log(0.01)/0.01; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05
hist(rexp(n=1000, rate = lambda_sigma))

P <- 3
mu_theta <- c(mean(log(obsY)),rep(0,P-1))
V_theta <- diag(c(1,rep(0.05,P-1)))

input <- list(N = nsize, K = nNeighbors, P = P, y = log(obsY), X = obsX, coords = obsCoords, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, mu_theta = mu_theta, V_theta = V_theta, a = ab[1], b = ab[2], lambda_sigma = lambda_sigma, lambda_tau = lambda_tau, sigma_multiplier = 0.05, tau_multiplier = 0.05)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/NNGP_HN.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          adapt_delta = 0.98,
                          max_treedepth = 12,
                          step_size = 0.25,
                          init = 1)
elapsed_time <- cmdstan_fit$time()
elapsed_time
elapsed_time$total/3600

cmdstan_fit$cmdstan_diagnose()
sampler_diag <- cmdstan_fit$sampler_diagnostics(format = "df")
str(sampler_diag)

## Posterior summaries
pars <- c("theta[1]","theta[2]","theta[3]","sigma","ell","tau")
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fit_summary %>% filter(variable %in% pars)
fixed_summary %>% print(digits = 3)

## Posterior draws
draws_df <- cmdstan_fit$draws(format = "df")
draws_df

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + facet_text(size = 15)

save(elapsed_time, fit_summary, fixed_summary, draws_df, file = paste0(fpath,"TemperatureDataAnalysis/logNNGP_Temps.RData"))

##################################################################
## Fitted value at each oberved sites
##################################################################
## Stan function exposed to be used 
source(paste0(fpath,"Rutilities/expose_cmdstanr_functions.R"))
exsf <- expose_cmdstanr_functions(model_path = stan_file)
args(exsf$vecchia_matern32_fitted_rng)


## Compute the means
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(obsX)
str(post_theta)
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
post_sigma <- as_tibble(draws_df) %>% .$sigma; str(post_sigma)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell <- as_tibble(draws_df) %>% .$ell; str(post_ell)

args(exsf$vecchia_matern32_fitted_rng)
l <- 1
logyfitted_list <- lapply(1:size_post_samples, function(l) exsf$vecchia_matern32_fitted_rng(y = log(obsY), mu = obsXtheta[l,], sigmasq = post_sigma[l]^2, tausq = post_tau[l]^2, lscale = post_ell[l], site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, neiID = lapply(1:nrow(neiMatInfo$NN_ind), function(i) neiMatInfo$NN_ind[i,]), N = nsize, K = nNeighbors))
yfitted_list <- lapply(logyfitted_list, exp)
yfitted_draws <- do.call(rbind,yfitted_list)
str(yfitted_draws)

yfitted_summary <- tibble(
  post.mean = apply(yfitted_draws, 2, mean),
  post.sd = apply(yfitted_draws, 2, sd),
  post.q2.5 = apply(yfitted_draws, 2, quantile2.5),
  post.q50 = apply(yfitted_draws, 2, quantile50),
  post.q97.5 = apply(yfitted_draws, 2, quantile97.5),
  y = obsY)
yfitted_summary

##################################################################
## Independent prediction at each predictions sites
##################################################################
## Stan function exposed to be used 
source(paste0(fpath,"Rutilities/expose_cmdstanr_functions.R"))
exsf <- expose_cmdstanr_functions(model_path = stan_file)
args(exsf$predict_nnnnglgc_rng)

size_post_samples <- nrow(draws_df); size_post_samples
## Compute the means
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(obsX)
str(post_theta)
l <- 1
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
prdXtheta <- t(sapply(1:size_post_samples, function(l) prdX %*% post_theta[l,])); str(prdXtheta)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_sigma <- as_tibble(draws_df) %>% .$sigma; str(post_sigma)
post_ell <- as_tibble(draws_df) %>% .$ell; str(post_ell)

## NNGP Preparation
psize <- nrow(prdCoords); psize
nei_info_pred <- FNN::get.knnx(obsCoords, prdCoords, k = nNeighbors); str(nei_info_pred)
pred2obsNeiID <- nei_info_pred$nn.index; str(pred2obsNeiID)
pred2obsDist <- nei_info_pred$nn.dist; str(pred2obsDist)

args(exsf$predict_vecchia_rng)
post_logypred <- exsf$predict_vecchia_rng(
  y = log(obsY), 
  obsX = obsX, 
  predX = prdX, 
  obsCoords = lapply(1:nrow(obsCoords), function(i) obsCoords[i,]),
  pred2obsDist = lapply(1:nrow(pred2obsDist), function(i) pred2obsDist[i,]), 
  pred2obsNeiID = lapply(1:nrow(pred2obsNeiID), function(i) pred2obsNeiID[i,]),
  beta = lapply(1:nrow(post_theta), function(i) post_theta[i,]), 
  sigma = post_sigma, 
  lscale = post_ell, 
  tau = post_tau, 
  nsize = nsize, 
  psize = psize, 
  postsize = size_post_samples)

post_ypred <- lapply(post_logypred, exp)
ypred_draws <- do.call(rbind,post_ypred); str(ypred_draws)
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
logs <- mean(logs_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,])); logs
CRPS <- mean(crps_sample(y = prdY[-id_full_missing], dat = t(ypred_draws)[-id_full_missing,])); CRPS

scores_df <- pred_summary %>% filter(!is.na(y)) %>%
  mutate(intervals = scoringutils::interval_score(true_values = y, lower = post.q2.5, upper = post.q50, interval_range = 0.95)) %>%
  mutate(btw = between(y,post.q2.5, post.q97.5)) %>%
  mutate(error = y - post.q50) %>%
  summarise(MAE = sqrt(mean(abs(error))), RMSE = sqrt(mean(error^2)), CVG = mean(btw),
            IS = mean(intervals)) %>%
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "logNNGP") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(elapsed_time, obsCoords, prdCoords, yfitted_summary, fit_summary, fixed_summary, draws_df, pred_summary, scores_df, file = paste0(fpath,"TemperatureDataAnalysis/logNNGP_Temps.RData"))

