matern32 <- function(d, sigma, lscale){
  ds <- sqrt(3)*d/lscale
  (sigma^2) * (1 + ds) * exp(-ds)
}
matern12 <- function(d, sigma, lscale){
  ds <- d/lscale
  (sigma^2) * exp(-ds)
}

############################################################################
## Prior elicitation for scale parameter
############################################################################
# P(lower_pars < pars < upper_pars) = probability

getIGamma <- function(par, lRange, uRange, prob) {
  
  lProb <- (1 - prob)/2
  uProb <- 1 - lProb
  
  invlRange <- 1/lRange
  invuRange <- 1/uRange
  
  c(qgamma(p = lProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invuRange, qgamma(p = uProb, shape = par[1], rate = par[2], lower.tail = TRUE) - invlRange)
}

############################################################################
## Frechet distribution
############################################################################

dfrechet <- function(x, alpha, sigma, log = FALSE){
  ldensity <- dweibull(1/x, shape = alpha, scale = 1/sigma, log = TRUE) - 2 * log(x)
  if (log){
    out <- ldensity
  } else{
    out <- exp(ldensity)
  }
  return(out)
}

qfrechet <- function(p, alpha, sigma, lower.tail = TRUE, log.p = FALSE){
  1/qweibull(1 - p, shape = alpha, scale = 1/sigma, lower.tail = lower.tail, log.p = log.p)
}

rfrechet <- function(n, alpha, sigma){
  1/rweibull(n, shape = alpha, scale = 1/sigma)
}

pfrechet <- function(q, alpha, sigma, lower.tail = TRUE, log.p = FALSE){
  pweibull(1/q, shape = alpha, scale = 1/sigma, lower.tail = !lower.tail, log.p = log.p)
}

############################################################################
## Inverse Gamma distribution
############################################################################
dinvgamma <- function (x, shape, scale, log = FALSE) {
  ldensity <- dgamma(1/x, shape, rate = scale, log = TRUE) - 2 * log(x)
  if (log){
    out <- ldensity
  } else{
    out <- exp(ldensity)
  }
  return(out)
}

qinvgamma <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE){
  1/qgamma(1 - p, shape, rate = scale, lower.tail = lower.tail, log.p = log.p)
}

rinvgamma <- function(n, shape, scale){
  1/rgamma(n, shape, rate = scale)
}

pinvgamma <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE){
  pgamma(1/q, shape, rate = scale, lower.tail = !lower.tail, log.p = log.p)
}


quantile2.5 <- function(x) quantile(x, prob = 0.025)
quantile50 <- function(x) quantile(x, prob = 0.50)
quantile97.5 <- function(x) quantile(x, prob = 0.975)


