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
source(paste0(fpath,"ExactVsApproximateMethods/gen_data_set4.R"))

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
# Preparing for Hilbert Space Approximate GP
################################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
quantile(obsDistVec, probs = c(1,2.5,52,50)/100)

ell_hat <- 0.5
c <- pmax(1.2, 4.5*ell_hat); c

m1 <- pmax(22,ceiling(3.42*c/ell_hat)); m1
m2 <- pmax(22,ceiling(3.42*c/ell_hat)); m2
mstar <- m1*m2; mstar

L <- c*Lstar; L
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
mu_theta <- c(mean(obsY),rep(0, P-1)); mu_theta
V_theta <- diag(c(5,rep(1,P-1))); V_theta
input <- list(N = nsize, M = mstar, P = P, y = obsY, X = obsX, coords = obsCoords, L = L, lambda = lambda, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, a = ab[1], b = ab[2], lambda_ell1 = lambda_ell1, lambda_ell2 = lambda_ell2, positive_skewness = 1, sigma1_multiplier = 0.50, sigma2_multiplier = 0.50, tau_multiplier = 0.25, gamma_multiplier = 0.40)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/HSHS_GLGC_HN.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 100,
                          iter_sampling = 100,
                          adapt_delta = 0.99,
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
pars <- c(paste0("theta[",1:P,"]"),"sigma1","sigma2","ell1","ell2","tau","gamma")
pars_true_df <- tibble(variable = pars, true = c(theta,sigma1,sigma2,lscale1,lscale2,tau,gamma))
fit_summary <- cmdstan_fit$summary(NULL, c("mean","sd","quantile50","quantile2.5","quantile97.5","rhat","ess_bulk","ess_tail"))
# cmdstan_fit$summary(variables = pars, posterior::default_summary_measures()[1:4], quantiles = ~ quantile(., probs = c(0.025, 0.975)), posterior::default_convergence_measures())
fixed_summary <- inner_join(pars_true_df, fit_summary)
fixed_summary %>% print(digits = 3)

## Posterior draws
draws_df <- cmdstan_fit$draws(format = "df")
class(draws_df)

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

## Recovery of z <- gamma*exp(z1) + z2
size_post_samples <- nrow(draws_df); size_post_samples
post_gamma <- as_tibble(draws_df) %>% .$gamma; str(post_gamma)
str(post_z1)
str(post_z2)
l <- 1
post_z <- t(sapply(1:size_post_samples, function(l) post_gamma[l]*exp(post_z1[l,]) + post_z2[l,])); str(post_z)
z <- gamma*exp(z1) + z2
z_summary <- tibble(z = z[idSampled],
                    post.mean = apply(post_z, 2, mean),
                    post.sd = apply(post_z, 2, sd),
                    post.q2.5 = apply(post_z, 2, quantile2.5),
                    post.q50 = apply(post_z, 2, quantile50),
                    post.q97.5 = apply(post_z, 2, quantile97.5))
z_summary
save(elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, z_summary, file = paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4.RData"))

#########################################################################
## Fitted value
#########################################################################

## Compute the means
size_post_samples <- nrow(draws_df); size_post_samples
post_theta <- as_tibble(draws_df) %>% select(starts_with("theta[")) %>% as.matrix() %>% unname(); str(post_theta)
str(obsX)
str(post_theta)
obsXtheta <- t(sapply(1:size_post_samples, function(l) obsX %*% post_theta[l,])); str(obsXtheta)
str(post_z1)
post_tau <- as_tibble(draws_df) %>% .$tau; str(post_tau)

l <- 1
yfitted_list <- lapply(1:size_post_samples, function(l){
  rnorm(n = nsize, mean = obsXtheta[l,] + post_gamma[l]*exp(post_z1[l,]) + post_z2[l,], sd = post_tau[l])
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

### Random effect z2 at predicted locations
psize <- nrow(prdCoords); psize
predH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda))); str(predH)
post_z2pred <- t(sapply(1:size_post_samples, function(l) predH %*% post_omega2[l,])); str(post_z2pred)

z2pred_summary <- tibble(
  z2 = prdZ2,
  post.mean = apply(post_z2pred, 2, mean),
  post.sd = apply(post_z2pred, 2, sd),
  post.q2.5 = apply(post_z2pred, 2, quantile2.5),
  post.q50 = apply(post_z2pred, 2, quantile50),
  post.q97.5 = apply(post_z2pred, 2, quantile97.5))
head(z2pred_summary)
mean(z2pred_summary[,"z2"] > z2pred_summary[,"post.q2.5"] & z2pred_summary[,"z2"] < z2pred_summary[,"post.q97.5"])

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
  mutate(ES = ES, logs = logs, CRPS = CRPS,  `Elapsed Time` = elapsed_time$total, Method = "HSHS_DataSet4") %>%
  select(Method,MAE,RMSE,CVG,CRPS,IS,ES,logs,`Elapsed Time`)
scores_df

save(sampler_diag, elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, z_summary, post_z, pred_summary, scores_df, file = paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4.RData"))

##################################################################
## Empirical density and CI
##################################################################
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
  geom_line(data = kde_df, aes(x = x, y = d, col = Key), linewidth = 0.35) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8),
        legend.title = element_blank())
ggsave(paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4_SpatialEffect_Density.png"), height = 4, width = 6)

save(sampler_diag, elapsed_time, fixed_summary, draws_df, kde_df, z_summary, post_z, yfitted_summary, pred_summary, scores_df, file = paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4.RData"))

##################################################################
## Empirical kernel density and Confidence Intervals
##################################################################

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
ggsave(paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4_SpatialEffect_Density.png"), height = 4, width = 6)

save(sampler_diag, kde_df, elapsed_time, fixed_summary, draws_df, z1_summary, z2_summary, z_summary, post_z, yfitted_summary, pred_summary, scores_df, file = paste0(fpath,"ExactVsApproximateMethods/HSHS_DataSet4.RData"))


