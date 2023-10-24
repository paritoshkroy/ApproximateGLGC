rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(nimble, warn.conflicts = FALSE)
library(coda)
library(nleqslv)

fname <- "GLGCMatern32NNNNAnalysis"

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
}
node_ind <- args[1]
# node_ind <- 201
rep_ind <- as.numeric(node_ind) - 200
rep_ind

##########################################################################
load("./CoordsData/coordsSimStudy.rda")
load("./GLGCMatern32Data/GLGCMatern32Data0.rda")
load(paste0("./GLGCMatern32Data/GLGCMatern32Data",rep_ind,".rda"))
cat("The",node_ind,"-th node is going to use",rep_ind,"-th replicated sample","\n")
##########################################################################

##########################################################################
source("Rfunctions/NNMatrix.R")
########################################################################
obsCoords <- coords[idSampled,]
str(obsCoords)
prdCoords <- coords[-idSampled,]
obsY <- y[idSampled]
prdY <- y[-idSampled]

obsDistMat <- fields::rdist(obsCoords)
str(obsDistMat)
obsDistVec <- obsDistMat[lower.tri(obsDistMat, diag = FALSE)]
obsMaxDist <- max(obsDistVec)
obsMedDist <- median(obsDistVec)
obsMinDist <- min(obsDistVec)

## NNGP preparation

nNeighbors <- 20
neiMatInfo <- NNMatrix(coords = obsCoords, n.neighbors = nNeighbors, n.omp.threads = 2)
str(neiMatInfo)
obsY <- obsY[neiMatInfo$ord] # ordered the data following neighborhood settings
obsCoords <- obsCoords[neiMatInfo$ord,] # ordered the data following neighborhood settings
obsZ1 <- z1[idSampled][neiMatInfo$ord]

# Constants
source("./Rfunctions/priorElicitationFunctions.R")
lLimit <- quantile(obsDistVec, prob = 0.025)
lLimit
uLimit <- quantile(obsDistVec, prob = 0.975)
uLimit

library(nleqslv)
hypIG <- nleqslv(c(3,1), getIGamma, lRange = lLimit, uRange = uLimit, prob = 0.98)
hypIG$x
plot(density(1/rgamma(n=10000, shape = hypIG$x[1], rate= hypIG$x[2])))

stan_data <- list(N = nsize, K = nNeighbors, y = obsY, neiID = neiMatInfo$NN_ind, site2neiDist = neiMatInfo$NN_dist, neiDistMat = neiMatInfo$NN_distM, hypShape = hypIG$x[1], hypScale = hypIG$x[2])
str(stan_data)

library(cmdstanr)
stan_file <- "GLGCMatern32NNNNb0.stan"
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

traceplot(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
stanfit_summary <- summary(stanfit, pars = c("beta","sigma1","lscale1","sigma2","lscale2","tau","gamma"))
summary_hyper <- stanfit_summary$summary    
summary_hyper

stanfit_summary <- summary(stanfit, pars = c("z"))
str(stanfit_summary$summary)
summary_z <- cbind(stanfit_summary$summary, z = obsZ1)
head(summary_z)
summary(unname(summary_z[,"Rhat"]))
mean(summary_z[,"z"] >= summary_z[,"2.5%"] & summary_z[,"z"] <= summary_z[,"97.5%"])

save(stanfit, elapsed_time, summary_hyper, summary_z, file = paste0("/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summary", fname, rep_ind,".rda"))

######################################################################
# Preparing for Prediction
######################################################################
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

post_beta <- list_of_post_samples$beta
attr(post_beta,"dimnames") <- NULL
post_beta <- t(post_beta)
str(post_beta)

post_z <- list_of_post_samples$z
attr(post_z,"dimnames") <- NULL
post_z <- t(post_z)
str(post_z)

n_post_samples <- length(post_tau)

obsX <- matrix(1, nrow = nsize, ncol = 1)
prdX <- matrix(1, nrow = nsite-nsize, ncol = 1)
str(obsX)
str(prdX)

Rcpp::sourceCpp("Cfunctions/predSumNNGLGCMatern32gamma.cpp")
pred_summary <- predSumNNGLGCMatern32gamma(obsY = stan_data$y, obsX = obsX, prdX = prdX, obsCoords = obsCoords, prdNeiDist = neiDist_pred, prdNeiID = neiID_pred, z = post_z, beta = post_beta, sigmasq1 = post_sigma1^2, lscale1 = post_lscale1^2, sigmasq2 = post_sigma2^2, lscale2 = post_lscale2, tausq = post_tau^2, gamma = post_gamma, iterprint = 500)
colnames(pred_summary) <- c("post.mean", "post.sd", "q2.5", "q50", "q97.5")
pred_summary <- cbind(pred_summary,y = prdY)
head(pred_summary)

save(elapsed_time, summary_hyper, summary_z, pred_summary, file = paste0("/home/pkroy/projects/def-aschmidt/pkroy/ApproximateGLGC/summary", fname, rep_ind,".rda"))



