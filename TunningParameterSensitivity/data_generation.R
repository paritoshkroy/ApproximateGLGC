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
#lscale1 <- sqrt(2)/2 
#lscale2 <- sqrt(2)/2
tau <- 0.5
tausq <- tau^2
gamma <- 1.5

distMat <- fields::rdist(coords)

SigmaX <- 1*0.25^abs(outer(1:2,1:2,'-'))
set.seed(10)
X <- cbind(1,cbind(rnorm(n=nsite),rnorm(n=nsite)) %*% t(chol(SigmaX)))
muX <- drop(X %*% theta)
z1 <- drop(crossprod(chol(matern32(d = fields::rdist(coords), sigma = sigma1, lscale = lscale1) + diag(x=1e-9, nrow = nsite, ncol = nsite)), rnorm(nsite)))
z2 <- drop(crossprod(chol(matern32(d = fields::rdist(coords), sigma = sigma2, lscale = lscale2) + diag(x=1e-9, nrow = nsite, ncol = nsite)), rnorm(nsite)))
linpred <- muX +  gamma * exp(z1) + z2
y <- rnorm(n = nsite, mean = linpred, sd = tau)
nsize <- 500
idSampled <- sample.int(n = nsite, size = nsize, replace = FALSE)
set.seed(NULL)
#hist(y, nclass = 21)
#hist(y[idSampled], nclass = 21)
#plot(density(y[idSampled]))
#boxplot(y[idSampled])

## Settings for HSGP following m = 3.42*c/(ell/S) with c>= 1.2
3.42*1.22/0.1
3.42*1.25/0.2
3.42*1.50/0.3
3.42*1.50/0.4
c1 <- c(1.22,1.25,rep(1.5,8)); c1
c2 <- c(1.22,1.25,rep(1.5,8)); c2
m1 <- c(42,23,rep(22,8)); m1
m2 <- c(42,23,rep(22,8)); m2

obj_all <- ls()
obj_keep <- c("idSampled", "y", "z1", "z2", "X", "theta", "sigma1", "sigma2", "lscale1", "lscale2", "tau", "coords", "gamma", "nsize", "m1", "m2", "c1", "c2")
obj_drop <- obj_all[!obj_all %in% c(obj_prev, obj_keep)]
rm(list = obj_drop)
gc()

