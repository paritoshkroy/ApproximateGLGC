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
nsize <- nrow(obsCoords); nsize
psize
obsY <- msst_df$temp; str(obsY)
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
## NNGP preparation
################################################################################
source(paste0(fpath,"Rutilities/NNMatrix.R"))
nNeighbors <- 5
neiMatInfo <- NNMatrix(coords = obsCoords, n.neighbors = nNeighbors, n.omp.threads = 2)
str(neiMatInfo)
obsY <- obsY[neiMatInfo$ord] # ordered the data following neighborhood settings
obsX <- obsX[neiMatInfo$ord,] # ordered the data following neighborhood settings
obsCoords <- obsCoords[neiMatInfo$ord,] # ordered the data following neighborhood settings

## Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit

lambda_sigma <- -log(0.01)/1; lambda_sigma
lambda_tau <- -log(0.01)/1; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05

library(nleqslv)
ab <- nleqslv(c(3,2), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), 0, uLimit)

P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
V_theta <- diag(c(10,rep(1,P-1)))

# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, K = nNeighbors, P = P, y = log(obsY), X = obsX, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma = lambda_sigma, lambda_tau = lambda_tau, a = ab[1], b = ab[2])
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/NNGP_Exp.stan")
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
pars <- c(paste0("theta[",1:P,"]"),"sigma","ell","tau")
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
fixed_summary <- fit_summary %>% filter(variable %in% pars)
fixed_summary %>% print(digits = 3)

## Posterior draws
draws_df <- cmdstan_fit$draws(format = "df")
draws_df

library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + facet_text(size = 15)

size_post_samples <- nrow(draws_df); size_post_samples
post_ell <- as_tibble(draws_df) %>% .$ell; str(post_ell)
post_sigma <- as_tibble(draws_df) %>% .$sigma; str(post_sigma)
save(elapsed_time, fixed_summary, draws_df, file = paste0(fpath,"SSTempDataAnalysis/logNNGP_GP_SST.RData"))

##################################################################
## Independent prediction at each predictions sites
##################################################################
source(paste0(fpath,"Rutilities/expose_cmdstanr_functions.R"))
exsf <- expose_cmdstanr_functions(model_path = stan_file)
args(exsf$predict_fullgp_rng)

## NNGP Preparation
psize <- nrow(prdCoords); psize
nei_info_pred <- FNN::get.knnx(obsCoords, prdCoords, k = nNeighbors); str(nei_info_pred)
pred2obsNeiID <- nei_info_pred$nn.index; str(pred2obsNeiID)
pred2obsDist <- nei_info_pred$nn.dist; str(pred2obsDist)

size_post_samples <- nrow(draws_df); size_post_samples
psize <- nrow(prdCoords); psize
post_sigma <- as_tibble(draws_df) %>% .$sigma; str(post_sigma)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell <- as_tibble(draws_df) %>% .$ell; str(post_ell)
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)

str(obsX)
str(post_theta)

obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
prdXtheta <- t(sapply(1:size_post_samples, function(l) prdX %*% post_theta[l,])); str(prdXtheta)
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

post_ypred <- lapply(post_logypred,exp)
ypred_draws <- do.call(rbind,post_ypred); str(ypred_draws)
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
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "logNNGP") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(elapsed_time, prdGrid, fixed_summary, draws_df, pred_summary, scores_df, file = paste0(fpath,"SSTempDataAnalysis/logNNGP_SST.RData"))

