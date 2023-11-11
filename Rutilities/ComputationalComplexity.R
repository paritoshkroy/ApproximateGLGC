## Computational Complexity
FullComputationalComplexity <- function(n) {
  2*n^3
}

NNNNComputationalComplexity <- function(n,m) {
  2*n*m^3
}

NNHSComputationalComplexity <- function(n,m,m1,m2) {
  n*m^3 + (n+1)*(m1*m2)
}

HSHSComputationalComplexity <- function(n,m1,m2) {
  2*(n+1)*(m1*m2) + n
}


NNNNComputationalComplexity(n = 50, m = 15) / HSHSComputationalComplexity(n = 50, m1 = 42, m2 = 42)

vector_lscale <- seq(0.05,0.65,l=13); vector_lscale
vector_c <- round(pmax(4.5*vector_lscale,1.2),2); vector_c
vector_m1 <- pmax(round(3.42*vector_c/vector_lscale,0),22); vector_m1
vector_m2 <- pmax(round(3.42*vector_c/vector_lscale,0),32); vector_m2
vector_m3 <- pmax(round(3.42*vector_c/vector_lscale,0),42); vector_m3
library(tidyverse)
setup <- data.frame(lscale = vector_lscale, c = vector_c, m1 = vector_m1, m2 = vector_m2, m3 = vector_m3)
setup <- setup %>%
  mutate(NN6_NN6 = NNNNComputationalComplexity(n = 500, m = 6)) %>%
  mutate(NN10_NN10 = NNNNComputationalComplexity(n = 500, m = 10)) %>%
  mutate(NN15_NN15 = NNNNComputationalComplexity(n = 500, m = 15)) %>%
  mutate(NN20_NN20 = NNNNComputationalComplexity(n = 500, m = 20)) %>%
  mutate(HSopt_HSopt1 = HSHSComputationalComplexity(n = 500, m1 = m1, m2 = m1)) %>%
  mutate(HSopt_HSopt2 = HSHSComputationalComplexity(n = 500, m1 = m2, m2 = m2)) %>%
  mutate(HSopt_HSopt3 = HSHSComputationalComplexity(n = 500, m1 = m3, m2 = m3))
setup %>% gather(Key,Value,-lscale,-c,-m1,-m2,-m3) %>%
  ggplot(aes(x = lscale, y = Value)) + 
  geom_point(col = "dimgray", shape = 1) +
  geom_path(aes(col = Key), linewidth = 0.7, alpha = 0.7) +
  theme_bw() +
  xlab(bquote("\u2113"[k])) +
  ylab("Computational Complexity") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        axis.title = element_text(size = 13))



lscale <- seq(0.05,1,l=10)
m <- 3.42*1.2/lscale;m
m1 <- c(80,36,30,25,22,22,22,22,22,22)
m2 <- c(80,36,30,25,22,22,22,22,22,22)
m <- c(15,15,15,15,15,10,10,10,10,10)
n <- 100000

3.42*1.5/0.1
sqrt(8)/2.75

RCC_list <- lapply(1:10, function(i) {
  FCC <- FullComputationalComplexity(n = n)
  NNNNCC <- NNNNComputationalComplexity(n = n, m = m[i])
  NNHSCC <- NNHSComputationalComplexity(n = n, m = m[i], m1= m1[i], m2 = m2[i])
  HSHSCC <- HSHSComputationalComplexity(n = n, m1 =m1[i], m2 = m2[i])
  tibble(Method = c("NNNN","NNHS","HSHS"), 
         `Relative Complexity` = c(NNNNCC/FCC, NNHSCC/FCC, HSHSCC/FCC),
         lscale = lscale[i])
})

RCC_df <- do.call(rbind, RCC_list) 
RCC_df

RCC_df %>% 
  ggplot(aes(x = lscale, y = `Relative Complexity`, group = Method)) + geom_line(aes(col = Method))


2/c(40,20,10,5,4,3,2)

m <- 3.42*1.5/seq(0.1,1,l=10)
m
3.42*1.2/0.08
m <- c(22,32,42,52)
ell <- seq(0.05,0.5,l=10)
sapply(1:10, function(i) m*ell[i]/3.42)
sapply(1:10, function(i) 4.5*ell[i])
sapply(1:10, function(i) pmax(pmin(m*ell[i]/3.42, 4.5*ell[i]),1.5))

ell <- seq(0.05,0.65,l=13); ell
c <- 4.5*ell; c
c <- round(pmax(4.5*ell,1.2),1); c
m <- round(3.42*c/ell,0); m
setup <- data.frame(lscale = ell, c = c, m = m)
setup

minimum_m <- function(c,l_by_S){
  ceiling(3.42*c/l_by_S)
}
m1.2 <- minimum_m(c=1.2,l_by_S = seq(0.1,1,l=10))
m1.3 <- minimum_m(c=1.3,l_by_S = seq(0.1,1,l=10))
m1.4 <- minimum_m(c=1.4,l_by_S = seq(0.1,1,l=10))
m1.5 <- minimum_m(c=1.5,l_by_S = seq(0.1,1,l=10))
data.frame(ell = seq(0.1,1,l=10), m1.2 = m1.2, m1.3 = m1.3, m1.4 = m1.4, m1.5 = m1.5)

plot(x = seq(0.1,1,l=10),y=m1.5, type= "b", col = "red")
lines(x = seq(0.1,1,l=10),y=m1.4, type= "b", col = "blue")
lines(x = seq(0.1,1,l=10),y=m1.3, type= "b", col = "green")
lines(x = seq(0.1,1,l=10),y=m1.2, type= "b", col = "magenta")

NNNN_GLGC_Exp 0.719 0.672 0.950 0.367 0.887  46.2 0.973         54395.
NNHS_GLGC_Exp 0.725 0.682 0.943 0.374 0.873  47.1 0.991         73260.
HSHS_GLGC_Exp 0.725 0.684 0.942 0.374 0.869  47.2 0.991          2007.