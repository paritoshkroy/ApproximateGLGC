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
node <- 30
fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"

##########################################################################
# ARC Preparation
##########################################################################
fpath <- "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/" #@ARC
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
node <- as.numeric(args[1])-250 ### specify correct node here
cat("The seed used to be ", node, "\n")

##########################################################################
# Data generation
##########################################################################
source(paste0(fpath,"Rutilities/utility_functions.R"))
source(paste0(fpath,"SmallLengthScale/data_generation.R"))

# partition as observed and predicted
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
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)
rm(obsDistMat)

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
obsZ1 <- z1[idSampled][neiMatInfo$ord]
obsZ2 <- z2[idSampled][neiMatInfo$ord]

################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat))); Lstar
c <- 1.20
m1 <- pmax(52,ceiling(3.42*c/(0.10/Lstar[1]))); m1
m2 <- pmax(52,ceiling(3.42*c/(0.10/Lstar[2]))); m2
mstar <- m1*m2; mstar
L <- c*Lstar; L
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)

## Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.99); uLimit

library(nleqslv)
ab <- nleqslv(c(5,0.1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
ab
curve(dinvgamma(x, shape = ab[1], scale = ab[2]), 0, uLimit)
summary(rinvgamma(n = 1000, shape = ab[1], scale = ab[2]))


## Exponential and PC prior
lambda_sigma1 <- -log(0.01)/1; lambda_sigma1
lambda_sigma2 <- -log(0.01)/1; lambda_sigma2
lambda_tau <- -log(0.01)/1; lambda_tau
pexp(q = 1, rate = lambda_tau, lower.tail = TRUE) ## P(tau > 1) = 0.05
lambda_ell1 <- as.numeric(-log(0.01)*lLimit); lambda_ell1
lambda_ell2 <- as.numeric(-log(0.01)*lLimit); lambda_ell2
pfrechet(q = lLimit, alpha = 1, sigma = lambda_ell2, lower.tail = TRUE) ## P(ell < lLimit) = 0.05
summary(rfrechet(n = 1000, alpha = 1, sigma = lambda_ell2))

## Stan input
head(obsX)
P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1)); mu_theta
V_theta <- diag(c(10,rep(1,P-1))); V_theta
# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, M = mstar, K = nNeighbors, P = P, y = obsY, X = obsX, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, coords = obsCoords, L = L, lambda = lambda, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, a = ab[1], b = ab[2], lambda_ell1 = lambda_ell1, lambda_ell2 = lambda_ell2, positive_skewness = 1, sigma1_multiplier = 0.50, sigma2_multiplier = 0.50, tau_multiplier = 0.25, gamma_multiplier = 0.50)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/NNHS_GLGC_HN.stan")
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
                          step_size = 0.25,
                          init = 1)
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

save(elapsed_time, fixed_summary, draws_df, z1_summary, post_z1, file = paste0(fpath,"SmallLengthScale/NNHS_GLGC",node,".RData"))

##################################################################
## Fitted value at each observed sites
##################################################################
## Stan function exposed to be used 
source(paste0(fpath,"Rutilities/expose_cmdstanr_functions.R"))
exsf <- expose_cmdstanr_functions(model_path = stan_file)
args(exsf$predict_nnhsglgc_rng)

args(exsf$vecchia_matern32_fitted_rng)
## Compute the means
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(obsX)
str(post_theta)
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
str(post_z1)
post_sigma2 <- as_tibble(draws_df) %>% .$sigma2; str(post_sigma2)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell2 <- as_tibble(draws_df) %>% .$ell2; str(post_ell2)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)

args(exsf$vecchia_matern32_fitted_rng)
l <- 1
yfitted_list <- lapply(1:size_post_samples, function(l){
  exsf$vecchia_matern32_fitted_rng(y = obsY, mu = obsXtheta[l,] + post_gamma[l]*exp(post_z1[l,]), sigmasq = post_sigma2[l]^2, tausq = post_tau[l]^2, lscale = post_ell2[l], site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, neiID = lapply(1:nrow(neiMatInfo$NN_ind), function(i) neiMatInfo$NN_ind[i,]), N = nsize, K = nNeighbors)
})

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

### Random effect z1 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z1pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega1[l,])); str(post_z1pred)

z1pred_summary <- tibble(
  z1 = prdZ1,
  post.mean = apply(post_z1pred, 2, mean),
  post.sd = apply(post_z1pred, 2, sd),
  post.q2.5 = apply(post_z1pred, 2, quantile2.5),
  post.q50 = apply(post_z1pred, 2, quantile50),
  post.q97.5 = apply(post_z1pred, 2, quantile97.5))
head(z1pred_summary)
mean(z1pred_summary[,"z1"] > z1pred_summary[,"post.q2.5"] & z1pred_summary[,"z1"] < z1pred_summary[,"post.q97.5"])

## Compute the means
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(post_z1)

str(obsX)
str(post_theta)

post_sigma1 <- as_tibble(draws_df) %>% .$sigma1; str(post_sigma1)
post_sigma2 <- as_tibble(draws_df) %>% .$sigma2; str(post_sigma2)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell1 <- as_tibble(draws_df) %>% .$ell1; str(post_ell1)
post_ell2 <- as_tibble(draws_df) %>% .$ell2; str(post_ell2)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(post_z1)

post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
str(post_z1pred)

## NNGP Preparation
psize <- nrow(prdCoords); psize
nei_info_pred <- FNN::get.knnx(obsCoords, prdCoords, k = nNeighbors); str(nei_info_pred)
pred2obsNeiID <- nei_info_pred$nn.index; str(pred2obsNeiID)
pred2obsDist <- nei_info_pred$nn.dist; str(pred2obsDist)

args(exsf$predict_nnhsglgc_rng)
post_ypred <- exsf$predict_nnhsglgc_rng(
  y = obsY, 
  obsX = obsX, 
  predX = prdX, 
  obsCoords = lapply(1:nrow(obsCoords), function(i) obsCoords[i,]),
  pred2obsDist = lapply(1:nrow(pred2obsDist), function(i) pred2obsDist[i,]), 
  pred2obsNeiID = lapply(1:nrow(pred2obsNeiID), function(i) pred2obsNeiID[i,]),
  beta = lapply(1:nrow(post_theta), function(i) post_theta[i,]), 
  z1 = lapply(1:nrow(post_z1), function(i) post_z1[i,]), 
  z1pred = lapply(1:nrow(post_z1pred), function(i) post_z1pred[i,]), 
  gamma = post_gamma, 
  sigma2 = post_sigma2, 
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
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "NNHS_GLGC") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(elapsed_time, fixed_summary, draws_df, z1_summary, post_z1, pred_summary, scores_df, file = paste0(fpath,"SmallLengthScale/NNHS_GLGC",node,".RData"))

##################################################################
## Recover latent vector z_2 at each observed sites
##################################################################
## Compute the means
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(obsX)
str(post_theta)
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
str(post_z1)
post_sigma2 <- as_tibble(draws_df) %>% .$sigma2; str(post_sigma2)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)
post_ell2 <- as_tibble(draws_df) %>% .$ell2; str(post_ell2)
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)

args(exsf$latent_matern32_rng)
post_z2_list <- lapply(1:size_post_samples, function(l){
  exsf$latent_matern32_rng(y = obsY, mu = obsXtheta[l,] + post_gamma[l]*exp(post_z1[l,]), sigma = post_sigma2[l], tau = post_tau[l], lscale = post_ell2[l], coords = lapply(1:nrow(obsCoords), function(i) obsCoords[i,]), N = nsize)
})

str(post_z2_list)
z2_draw <- do.call(rbind, post_z2_list)
str(z2_draw)
z2_summary <- tibble(
  post.mean = apply(z2_draw, 2, mean),
  post.sd = apply(z2_draw, 2, sd),
  post.q2.5 = apply(z2_draw, 2, quantile2.5),
  post.q50 = apply(z2_draw, 2, quantile50),
  post.q97.5 = apply(z2_draw, 2, quantile97.5),
  z2 = obsZ2)
z2_summary %>% summarise(CVG =  mean(between(x = z2, left = post.q2.5, right = post.q97.5)))
z2_summary %>% ggplot(aes(x = post.q50, y = z2)) + geom_point()

## Obtain z
z_draw_list <- lapply(1:size_post_samples, function(l) post_gamma[l]*exp(post_z1[l,]) + z2_draw[l,])
z_draw <- do.call(rbind, z_draw_list)
z_summary <- tibble(
  post.mean = apply(z_draw, 2, mean),
  post.sd = apply(z_draw, 2, sd),
  post.q2.5 = apply(z_draw, 2, quantile2.5),
  post.q50 = apply(z_draw, 2, quantile50),
  post.q97.5 = apply(z_draw, 2, quantile97.5),
  z = gamma*exp(obsZ1) + obsZ2)
z_summary %>% summarise(CVG =  mean(between(x = z, left = post.q2.5, right = post.q97.5)))
z_summary %>% ggplot(aes(x = post.q50, y = z)) + geom_point()

ggplot(z_summary) + 
  geom_density(aes(x = post.mean, col = "Posterior mean")) + 
  geom_density(aes(x = z, col = "True value")) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8),
        legend.title = element_blank())

ggplot(z_summary, aes(x = 1:nrow(z_summary))) + 
  geom_point(aes(y = post.mean)) + 
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25) +
  geom_point(aes(y = z), col = "red") +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank())

# R(K) for a normal
Rk <- 1 / (2 * sqrt(pi))

# Compute the kde (NR bandwidth)
kdeObs <- density(z_summary$z, from = min(z_summary$z), to = max(z_summary$z), n = nsize, bw = "nrd")
kdePm <- density(z_summary$post.mean, from = min(z_summary$post.mean), to = max(z_summary$post.mean), n = nsize, bw = "nrd")

# Selected bandwidth
hObs <- kdeObs$bw
hPm <- kdePm$bw

# Estimate the variance
var_kdeObs_hat <- kdeObs$y * Rk / (nsize * hObs)
var_kdePm_hat <- kdePm$y * Rk / (nsize * hPm)

# CI with estimated variance
alpha <- 0.05
z_alpha2 <- qnorm(1 - alpha / 2)
lciObs <- kdeObs$y - z_alpha2 * sqrt(var_kdeObs_hat)
uciObs <- kdeObs$y + z_alpha2 * sqrt(var_kdeObs_hat)

lciPm <- kdePm$y - z_alpha2 * sqrt(var_kdePm_hat)
uciPm <- kdePm$y + z_alpha2 * sqrt(var_kdePm_hat)

# Plot estimate, CIs and expectation
kdeObs_df <- tibble(x = kdeObs$x, d = kdeObs$y, lci = lciObs, uci = uciObs) %>% mutate(Key = 1)
kdePm_df <- tibble(x = kdePm$x, d = kdePm$y, lci = lciPm, uci = uciPm) %>% mutate(Key = 2)
kde_df <- rbind(kdeObs_df,kdePm_df)
kde_df <- kde_df %>% mutate(Key = factor(Key, labels = c("True","Estimated")))
ggplot(z_summary) +
  geom_histogram(aes(x = z, y = after_stat(density)), 
                 bins = 21, fill = NA, col = "dimgray") +
  geom_ribbon(data = kde_df, aes(x = x, ymin = lci, ymax = uci, fill = Key), alpha = 0.5) +
  geom_line(data = kde_df, aes(x = x, y = d, col = Key), 
            linetype = "dashed", linewidth = 0.5) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8),
        legend.title = element_blank())
ggsave(paste0(fpath,"SmallLengthScale/NNHS_GLGC_SPDensity",node,".png"), height = 4, width = 6)

save(sampler_diag, elapsed_time, fixed_summary, draws_df, pred_summary, scores_df, yfitted_summary, z_summary, kde_df, file = paste0(fpath,"SmallLengthScale/NNHS_GLGC",node,".RData"))

