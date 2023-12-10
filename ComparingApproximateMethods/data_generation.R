obj_prev <- ls()

library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

coords <- unname(as.matrix(expand.grid(x = seq(-0.99, 0.99, length.out = 100), y = seq(-0.99, 0.99, length.out = 100))))
nsite <- nrow(coords)
theta <- c(5,-1,1)
sigma1 <- 1
sigma2 <- 1
sigma1sq <- sigma1^2
sigma2sq <- sigma2^2
lscale1 <- 0.35
lscale2 <- 0.35
tau <- 0.5
tausq <- tau^2
gamma <- 1.75

distMat <- fields::rdist(coords)

SigmaX <- 1*0.25^abs(outer(1:2,1:2,'-'))
set.seed(350) # seed for generating Covariates
X <- cbind(1,cbind(rnorm(n=nsite),rnorm(n=nsite)) %*% t(chol(SigmaX)))
muX <- drop(X %*% theta)
set.seed(NULL)

set.seed(node*350) # seed for generating the random effect and response
z1 <- drop(crossprod(chol(matern32(d = fields::rdist(coords), sigma = sigma1, lscale = lscale1) + diag(x=1e-9, nrow = nsite, ncol = nsite)), rnorm(nsite)))
z2 <- drop(crossprod(chol(matern32(d = fields::rdist(coords), sigma = sigma2, lscale = lscale2) + diag(x=1e-9, nrow = nsite, ncol = nsite)), rnorm(nsite)))
linpred <- muX +  gamma * exp(z1) + z2
y <- rnorm(n = nsite, mean = linpred, sd = tau)
set.seed(NULL)

set.seed(1000) # seed for selecting sampled locations
nsize <- 1000
idSampled <- sample.int(n = nsite, size = nsize, replace = FALSE)
set.seed(NULL)

obj_all <- ls()
obj_keep <- c("idSampled", "y", "z1", "z2", "X", "theta", "sigma1", "sigma2", "lscale1", "lscale2", "tau", "coords", "gamma", "nsize")
obj_drop <- obj_all[!obj_all %in% c(obj_prev, obj_keep)]
rm(list = obj_drop)
gc()

