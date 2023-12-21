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

## Prior elicitation
lLimit <- quantile(obsDistVec, prob = 0.01); lLimit
uLimit <- quantile(obsDistVec, prob = 0.50); uLimit

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
P <- 3
mu_theta <- c(mean(obsY),rep(0,P-1))
V_theta <- diag(c(10,rep(1,P-1)))

# Keep in mind that the data should be ordered following nearest neighbor settings
input <- list(N = nsize, K = nNeighbors, P = P, y = obsY, X = obsX, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, mu_theta = mu_theta, V_theta = V_theta, lambda_sigma1 = lambda_sigma1, lambda_sigma2 = lambda_sigma2, lambda_tau = lambda_tau, a = ab[1], b = ab[2], lambda_ell1 = lambda_ell1, lambda_ell2 = lambda_ell2, positive_skewness = 0, sigma1_multiplier = 1, sigma2_multiplier = 1, tau_multiplier = 1)
str(input)

library(cmdstanr)
stan_file <- paste0(fpath,"StanFiles/NNNN_GLGC_HN.stan")
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
elapsed_time # 0.516789
cmdstan_fit$output()
# "Gradient evaluation took 0.031852 seconds"
# "1000 transitions using 10 leapfrog steps per transition would take 318.52 seconds."                
