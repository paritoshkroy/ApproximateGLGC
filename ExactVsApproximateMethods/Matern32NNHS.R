rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

fname <- "Matern32NNHS"
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
## NNGP preparation
################################################################################
source("Rfunctions/NNMatrix.R")
nNeighbors <- 20
neiMatInfo <- NNMatrix(coords = obsCoords, n.neighbors = nNeighbors, n.omp.threads = 2)
str(neiMatInfo)
obsY <- obsY[neiMatInfo$ord] # ordered the data following neighborhood settings
obsX <- obsX[neiMatInfo$ord,] # ordered the data following neighborhood settings
obsCoords <- obsCoords[neiMatInfo$ord,] # ordered the data following neighborhood settings
obsZ1 <- z1[idSampled][neiMatInfo$ord] # ordered the data following neighborhood settings

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

# Constants
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
V_beta <- diag(c(2.5*sd(obsY),rep(1,P)))
# Keep in mind that the data should be ordered following nearest neighbor settings
stan_data <- list(N = nsize, M = mstar, K = nNeighbors, P = 2, y = obsY, X = obsX, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, coords = obsCoords, L = L, lambda = lambda, mu_beta = mu_beta, V_beta = V_beta, hypShape = hypIG$x[1], hypScale = hypIG$x[2])
str(stan_data)

library(cmdstanr)
stan_file <- "GLGCMatern32NNHSXb.stan"
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
max(rowSums(elapsed_time))/60

save(stanfit, elapsed_time, file = "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summaryMatern32NNHS.rda")

sampler_params <- get_sampler_params(stanfit, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)

rstan::traceplot(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
rstan::traceplot(stanfit, pars = paste0("omega[",sample.int(n = mstar, size = 6, replace = FALSE),"]"))

stanfit_summary <- summary(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
summary_hyper <- stanfit_summary$summary    
summary_hyper

stanfit_summary <- summary(stanfit, pars = c("omega"))
str(stanfit_summary$summary)
head(stanfit_summary$summary)
summary(unname(stanfit_summary$summary[,"Rhat"]))

### Recovering the latent process z1
omega_samples <- rstan::extract(stanfit, pars = "omega")$omega
str(omega_samples)
eigenfunction_compute <- function(x, L, lambda) { 
  apply(sqrt(1/L) * sin(sqrt(lambda) %*% diag(x + L)), 1, prod)
}
obsH <- t(apply(obsCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(obsH)
z1_samples <- obsH %*% t(omega_samples)
str(z1_samples)
summary_z1 <- cbind(
  post.mean = apply(z1_samples, 1, mean),
  post.sd = apply(z1_samples, 1, sd),
  q2.5 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.025)),
  q50 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.5)),
  q97.5 = apply(z1_samples, 1, function(x) quantile(x, prob = 0.975)),
  z = obsZ1)
head(summary_z1)
mean(summary_z1[,"z"] > summary_z1[,"q2.5"] & summary_z1[,"z"] < summary_z1[,"q97.5"])

data.frame(summary_z1) %>%
  ggplot(aes(x = 1:nsize)) + 
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), linewidth = 0.25) +
  geom_point(aes(y = z), col = 2, size = 0.7, shape = 1)


### Calculation for predicting responses at predicted locations
psize <- nrow(prdCoords)
psize
prdH <- t(apply(prdCoords, 1, function(x) eigenfunction_compute(x, L = L, lambda = lambda)))
str(prdH)

nei_info_pred <- FNN::get.knnx(obsCoords, prdCoords, k = nNeighbors)
str(nei_info_pred)
neiID_pred <- nei_info_pred$nn.index
str(neiID_pred)
neiDist_pred <- nei_info_pred$nn.dist
str(neiDist_pred)

list_of_post_samples <- extract(stanfit)
print(names(list_of_post_samples))

post_sigma1 <- list_of_post_samples$sigma1
attr(post_sigma1,"dimnames") <- NULL
str(post_sigma1)

post_sigma2 <- list_of_post_samples$sigma2
attr(post_sigma2,"dimnames") <- NULL
str(post_sigma2)

post_lscale1 <- list_of_post_samples$lscale1
attr(post_lscale1,"dimnames") <- NULL
str(post_lscale1)

post_lscale2 <- list_of_post_samples$lscale2
attr(post_lscale2,"dimnames") <- NULL
str(post_lscale2)

post_tau <- list_of_post_samples$tau
attr(post_tau,"dimnames") <- NULL
str(post_tau)

post_gamma <- list_of_post_samples$gamma
attr(post_gamma,"dimnames") <- NULL
str(post_gamma)

post_omega <- list_of_post_samples$omega
attr(post_omega,"dimnames") <- NULL
post_omega <- t(post_omega)
str(post_omega)

post_beta <- list_of_post_samples$beta
attr(post_beta,"dimnames") <- NULL
post_beta <- t(post_beta)
str(post_beta)

n_post_samples <- length(post_tau)

str(post_omega)
str(obsH)
str(obsX)
str(prdX)
str(post_beta)
str(obsH)
str(post_omega)
post_linpred <- obsX %*% post_beta + sapply(1:length(post_gamma), function(i) post_gamma[i]*exp(obsH %*% post_omega[,i]))
str(post_linpred)

pred_linpred <- prdX %*% post_beta + sapply(1:length(post_gamma), function(i) post_gamma[i]*exp(prdH %*% post_omega[,i]))
str(pred_linpred)

## given linear predicted value X * beta + gamma*exp(z) prediction is the same as the marginal NNGP
Rcpp::sourceCpp("Cfunctions/predSumNNGPMatern32.cpp")
pred_summary <- predSumNNGPMatern32(obsY = stan_data$y, obsMuY = post_linpred, prdMuY = pred_linpred, obsCoords  = obsCoords, prdNeiDist = neiDist_pred, prdNeiID = neiID_pred, sigmasq = post_sigma2^2, lscale = post_lscale2, tausq = post_tau^2, iterprint = 500) 
colnames(pred_summary) <- c("post.mean", "post.sd", "q2.5", "q50", "q97.5")
pred_summary <- cbind(pred_summary, y = prdY)
head(pred_summary)
save(elapsed_time, summary_hyper, summary_z1, pred_summary, file = "/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summaryMatern32NNHS.rda")
