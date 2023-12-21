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

######################################################################
# Generating the data
######################################################################
source(paste0(fpath,"Rutilities/utility_functions.R"))
source(paste0(fpath,"GradientEvaluatingTime/gen_data_set1.R"))

######################################################################
# partition as observed and predicted
######################################################################
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
obsMaxDist <- max(obsDistVec); obsMaxDist
obsMedDist <- median(obsDistVec); obsMedDist
obsMinDist <- min(obsDistVec); obsMinDist
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
obsZ1 <- obsZ1[neiMatInfo$ord]
obsZ2 <- obsZ2[neiMatInfo$ord]

################################################################################
# Preparing for Hilbert Space Approximate GP
################################################################################
xRangeDat <- c(-1,1)
yRangeDat <- c(-1,1)
Lstar <- c(max(abs(xRangeDat)), max(abs(yRangeDat)))
quantile(obsDistVec, probs = c(1,2.5,52,50)/100)

ell_hat <- 0.2
c <- 1 + 2*ell_hat; c
c <- pmax(1.2, 4.5*ell_hat); c

m1 <- pmax(32,ceiling(3.42*c/ell_hat)); m1
m2 <- pmax(32,ceiling(3.42*c/ell_hat)); m2
mstar <- m1*m2; mstar

L <- c*Lstar; L
S <- unname(as.matrix(expand.grid(S2 = 1:m1, S1 = 1:m2)[,2:1]))
str(S)
lambda <- ((pi*S)/(2*L))^2
str(lambda)
head(lambda)

## Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit

library(nleqslv)
ab <- nleqslv(c(5,1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)$x
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
P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
V_theta <- diag(c(10,rep(1,P-1)))

# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, M = mstar, K = nNeighbors, P = P, y = obsY, X = obsX, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, coords = obsCoords, L = L, lambda = lambda, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, a = ab[1], b = ab[2], lambda_ell1 = lambda_ell1, lambda_ell2 = lambda_ell2, positive_skewness = 0, sigma1_multiplier = 1, sigma2_multiplier = 1, tau_multiplier = 1)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/NNHS_GLGC_HN.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)
mod$check_syntax(pedantic = TRUE)
mod$print()
cmdstan_fit <- mod$sample(data = input, 
                          chains = 1,
                          parallel_chains = 1,
                          iter_warmup = 1,
                          iter_sampling = 1,
                          adapt_delta = 0.99,
                          max_treedepth = 15,
                          step_size = 0.25)
elapsed_time <- cmdstan_fit$time()
elapsed_time # 1.532892
cmdstan_fit$output()
# "Gradient evaluation took 0.017859 seconds" 
# "1000 transitions using 10 leapfrog steps per transition would take 178.59 seconds."

