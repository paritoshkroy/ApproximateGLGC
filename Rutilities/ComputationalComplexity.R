rm(list =ls())
library(tidyverse)
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

## Combined NNNN NNHS and HSHS
## c = 1.20
vector_lscale <- seq(0.05,0.5,l=10); vector_lscale
round(4.5*vector_lscale,2)
vector_c <- pmax(1.20, round(4.5*vector_lscale,2)); vector_c
#vector_c <- 1 + 1.5*vector_lscale; vector_c
cbind(vector_lscale, vector_c, m = ceiling(3.42*vector_c/vector_lscale))
vector_m1 <- pmax(16,ceiling(3.42*vector_c/vector_lscale)); vector_m1
vector_m2 <- pmax(17,ceiling(3.42*vector_c/vector_lscale)); vector_m2
vector_m3 <- pmax(21,ceiling(3.42*vector_c/vector_lscale)); vector_m3
vector_m4 <- pmax(28,ceiling(3.42*vector_c/vector_lscale)); vector_m4
vector_m5 <- pmax(42,ceiling(3.42*vector_c/vector_lscale)); vector_m5
vector_m6 <- pmax(83,ceiling(3.42*vector_c/vector_lscale)); vector_m6


library(tidyverse)
setup <- tibble(lscale = vector_lscale, c = vector_c, m1 = vector_m1, m2 = vector_m2, m3 = vector_m3, m4 = vector_m4, m5 = vector_m5, m6 = vector_m6)
setup <- setup %>% gather(x,m,-c,-lscale) %>% select(-x) %>% distinct()
setup
table(setup$m)
setup %>% filter(lscale == 0.1)

setup <- setup %>% mutate(HSHS = HSHSComputationalComplexity(n = 1000, m1 = m, m2 = m)) %>%
  mutate(NNHS5 = NNHSComputationalComplexity(n = 1000, m = 5, m1 = m, m2 = m)) %>%
  mutate(NNHS10 = NNHSComputationalComplexity(n = 1000, m = 10, m1 = m, m2 = m)) %>%
  mutate(NNHS15 = NNHSComputationalComplexity(n = 1000, m = 15, m1 = m, m2 = m)) %>%
  mutate(NNHS20 = NNHSComputationalComplexity(n = 1000, m = 20, m1 = m, m2 = m))
setup
setup <- setup %>% mutate(Method = "HSHS")
setup
setup %>% filter(lscale == 0.50)
setup %>%
  filter(lscale %in% c(0.05,0.10,0.20,0.30,0.40,0.50)) %>%
  ggplot(aes(x = m, y = HSHS)) + 
  geom_path(linewidth = 0.7, color = "blue") +
  geom_path(aes(y = NNHS5), linetype = "dotted", color = "green", linewidth = 0.7) +
  geom_path(aes(y = NNHS10), linetype = "dotdash", color = "green",  linewidth = 0.7) +
  geom_path(aes(y = NNHS15), linetype = "dashed", color = "green",  linewidth = 0.7) +
  geom_path(aes(y = NNHS20), linetype = "longdash", color = "green", linewidth = 0.7) +
  geom_point(size = 1, shape = 20) +
  geom_point(aes(y = NNHS5), size =1, shape = 1) +
  geom_point(aes(y = NNHS10), size =1, shape = 2) +
  geom_point(aes(y = NNHS15), size =1, shape = 3) +
  geom_point(aes(y = NNHS20), size =1, shape = 4) +
  facet_wrap(~lscale, nrow = 2, labeller = label_bquote("\u2113"[k]~"="~.(lscale))) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 5), linetype = "1"), color = "red", linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 10), linetype = "2"), color = "red", linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 15), linetype = "3"), color = "red", linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 20), linetype = "4"), color = "red", linewidth = 0.5) +
  scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash"),
                        labels = c("m = 5","m = 10","m = 15","m = 20")) +
  theme_bw() +
  xlab(bquote(Number~of~basis~functions~under~HSHS)) +
  ylab(bquote(Computational~Cost)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.07,0.85),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank())
ggsave(filename = "./Rutilities/Computational_Complexity_Plots.png", height = 6, width = 9)

