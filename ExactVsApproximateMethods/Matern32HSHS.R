rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

fname <- "Matern32HSHS"
load("./DataCompared2FullGLGC/GLGCDataMatern32.rda")

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

################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
xRangeDat <- sDomain[,1]
yRangeDat <- sDomain[,2]
m1 <- 22; m2 <- 22; mstar <- m1*m2
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
c <- c(1.5,1.5)
L <- c*Lstar
str(L)
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)


source("./Rfunctions/priorElicitationFunctions.R")
lLimit <- quantile(obsDistVec, prob = 0.025)
lLimit
uLimit <- quantile(obsDistVec, prob = 0.975)
uLimit

library(nleqslv)
hypIG <- nleqslv(c(3,1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)
hypIG$x
plot(density(1/rgamma(n=1000, shape = hypIG$x[1], rate = hypIG$x[2])))

head(obsX)
P <- 2
mu_beta <- c(mean(obsY),rep(0, P))
V_beta <- diag(c(2.5*var(obsY),rep(1,P)))
stan_data <- list(N = nsize, M = mstar, P = 2, y = obsY, X = obsX, coords = obsCoords, L = L, lambda = lambda, mu_beta = mu_beta, V_beta = V_beta, hypShape = hypIG$x[1], hypScale = hypIG$x[2])
str(stan_data)

library(cmdstanr)
stan_file <- "GLGCMatern32HSHSXb.stan"
mod <- cmdstan_model(stan_file)
mod$print()
cmdstan_fit <- mod$sample(data = stan_data, 
                          chains = 4,
                          parallel_chains = 4,
                          iter_warmup = 1500,
                          iter_sampling = 500,
                          adapt_delta = 0.99,
                          max_treedepth = 12,
                          step_size = 0.5)
cmdstan_fit$time()
cmdstan_fit$cmdstan_diagnose()

library(rstan)
stanfit <- read_stan_csv(cmdstan_fit$output_files())

elapsed_time <- get_elapsed_time(stanfit)
print(elapsed_time)

save(stanfit, elapsed_time, file = "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summaryMatern32HSHS.rda")

load("~/McGill/ApproximateGLGC/summaryMatern32HSHS.rda")

sampler_params <- get_sampler_params(stanfit, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)

rstan::traceplot(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
rstan::traceplot(stanfit, pars = paste0("omega1[",sample.int(n = mstar, size = 6, replace = FALSE),"]"))
rstan::traceplot(stanfit, pars = paste0("omega2[",sample.int(n = mstar, size = 6, replace = FALSE),"]"))

stanfit_summary <- summary(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
summary_hyper <- stanfit_summary$summary    
summary_hyper

stanfit_summary <- summary(stanfit, pars = c("omega1"))
str(stanfit_summary$summary)
head(stanfit_summary$summary)
summary(unname(stanfit_summary$summary[,"Rhat"]))

stanfit_summary <- summary(stanfit, pars = c("omega2"))
str(stanfit_summary$summary)
head(stanfit_summary$summary)
summary(unname(stanfit_summary$summary[,"Rhat"]))

### Recovering the latent process z1
omega1_samples <- rstan::extract(stanfit, pars = "omega1")$omega1
str(omega1_samples)
eigenfunction_compute <- function(x, L, lambda) { 
  apply(sqrt(1/L) * sin(sqrt(lambda) %*% diag(x + L)), 1, prod)
}
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
z1_samples <- obsH %*% t(omega1_samples)
str(z1_samples)
summary_z1 <- cbind(
  post.mean = apply(z1_samples, 1, mean),
  post.sd = apply(z1_samples, 1, sd),
  q2.5 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.5)),
  q97.5 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.975)),
  z1 = z1[idSampled])
head(summary_z1)
mean(summary_z1[,"z1"] > summary_z1[,"q2.5"] & summary_z1[,"z1"] < summary_z1[,"q97.5"])

data.frame(summary_z1) %>%
  ggplot(aes(x = 1:nsize)) + 
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.25) +
  geom_point(aes(y = z1), col = 2, size = 0.7, shape = 1)


## For z2
omega2_samples <- rstan::extract(stanfit, pars = "omega2")$omega2
str(omega2_samples)
z2_samples <- obsH %*% t(omega2_samples)
str(z2_samples)
summary_z2 <- cbind(
  post.mean = apply(z2_samples, 1, mean),
  post.sd = apply(z2_samples, 1, sd),
  q2.5 = apply(z2_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(z2_samples, 1, function(x) quantile(x, prob = 0.50)),
  q97.5 = apply(z2_samples, 1, function(x) quantile(x, prob = 0.975)),
  z2 = z2[idSampled])
head(summary_z2)
mean(summary_z2[,"z2"] > summary_z2[,"q2.5"] & summary_z2[,"z2"] < summary_z2[,"q97.5"])

data.frame(summary_z2) %>%
  ggplot(aes(x = 1:nsize)) + 
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.25) +
  geom_point(aes(y = z2), col = 2, size = 0.7, shape = 1)


### Calculation for predicting responses at predicted locations
psize <- nrow(prdCoords)
psize
prdH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(prdH)

prd_z1_samples <- prdH %*% t(omega1_samples)
str(prd_z1_samples)
pred_summary_z1 <- cbind(
  post.mean = apply(prd_z1_samples, 1, mean),
  post.sd = apply(prd_z1_samples, 1, sd),
  q2.5 = apply(prd_z1_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(prd_z1_samples, 1, function(x) quantile(x, prob = 0.5)),
  q97.5 = apply(prd_z1_samples, 1, function(x) quantile(x, prob = 0.975)),
  z1 = z1[-idSampled])
head(pred_summary_z1)
mean(pred_summary_z1[,"z1"] > pred_summary_z1[,"q2.5"] & pred_summary_z1[,"z1"] < pred_summary_z1[,"q97.5"])

data.frame(pred_summary_z1) %>%
  ggplot(aes(x = 1:psize)) + 
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.25) +
  geom_point(aes(y = z1), col = 2, size = 0.7, shape = 1)

## For z2
prd_z2_samples <- prdH %*% t(omega2_samples)
str(prd_z2_samples)
pred_summary_z2 <- cbind(
  post.mean = apply(prd_z2_samples, 1, mean),
  post.sd = apply(prd_z2_samples, 1, sd),
  q2.5 = apply(prd_z2_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(prd_z2_samples, 1, function(x) quantile(x, prob = 0.50)),
  q97.5 = apply(prd_z2_samples, 1, function(x) quantile(x, prob = 0.975)),
  z2 = z2[-idSampled])
head(pred_summary_z2)
mean(pred_summary_z2[,"z2"] > pred_summary_z2[,"q2.5"] & pred_summary_z2[,"z2"] < pred_summary_z2[,"q97.5"])

data.frame(pred_summary_z2) %>%
  ggplot(aes(x = 1:psize)) + 
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.25) +
  geom_point(aes(y = z2), col = 2, size = 0.7, shape = 1)

post_beta_samples <- t(rstan::extract(stanfit, pars = "beta")$beta)
attr(post_beta_samples, "dimnames") <- NULL
str(post_beta_samples)
str(prdX)
prd_Xbeta_samples <- prdX %*% post_beta_samples
str(prd_Xbeta_samples)
str(prd_z1_samples)
str(prd_z2_samples)
post_tau_samples <- rstan::extract(stanfit, pars = "tau")$tau
attr(post_tau_samples, "dimnames") <- NULL
str(post_tau_samples)
post_gamma_samples <- rstan::extract(stanfit, pars = "gamma")$gamma
attr(post_gamma_samples, "dimnames") <- NULL
str(post_gamma_samples)

npost_samples <- length(post_tau_samples)
npost_samples
prd_y_samples <- sapply(1:npost_samples, function(i) prd_Xbeta_samples[,i] + post_gamma_samples[i] * exp(prd_z1_samples[,i]) + prd_z2_samples[,i] + rnorm(n = psize, mean = 0, sd = post_tau_samples[i]))
str(prd_y_samples)

pred_summary <- cbind(
  post.mean = apply(prd_y_samples, 1, mean),
  post.sd = apply(prd_y_samples, 1, sd),
  q2.5 = apply(prd_y_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(prd_y_samples, 1, function(x) quantile(x, prob = 0.50)),
  q97.5 = apply(prd_y_samples, 1, function(x) quantile(x, prob = 0.975)),
  y = y[-idSampled])
head(pred_summary)
mean(pred_summary[,"y"] > pred_summary[,"q2.5"] & pred_summary[,"y"] < pred_summary[,"q97.5"])

as.data.frame(pred_summary) %>%
  ggplot(aes(x= post.mean, y = y)) + geom_point(size = 0.24) +
  tune::coord_obs_pred() +
  xlab("Observed Response") +
  ylab("Mean of Posterior Predictive Distribution") +
  theme_bw()

save(elapsed_time, summary_hyper, summary_z1, summary_z2, pred_summary, pred_summary_z1, pred_summary_z2, file = "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summaryMatern32HSHS.rda")

