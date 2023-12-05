dglgc <- function(z, sigma1, sigma2, gamma) {
  fdensity <- integrate(f = function(x) {
      (1/(2*pi*sigma1*sigma2)*x) * 
        exp(- log(x)^2/(2*sigma1^2) - (z - gamma * x)^2/(2*sigma2^2))
    }, 0, Inf)$value
  
  return(fdensity)
}

library(ggplot2)
ggplot() + xlim(-19,2) + 
  geom_function(aes(colour = "0.25", linetype = "0.75"), 
                linewidth = 0.25, fun = Vectorize(dglgc),
                args = list(sigma1 = 0.90, sigma2 = 1, gamma = -1.03))

  