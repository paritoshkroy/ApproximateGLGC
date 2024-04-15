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
vector_lscale <- c(0.05,0.10,0.20,0.30,0.40,0.50); vector_lscale
round(4.5*vector_lscale,2)
vector_c <- pmax(1.20, round(4.5*vector_lscale,2)); vector_c
cbind(vector_lscale, vector_c, m = pmax(22,ceiling(3.42*vector_c/vector_lscale)))
vector_m1 <- pmax(22,ceiling(3.42*vector_c/vector_lscale)); vector_m1
vector_m2 <- pmax(32,ceiling(3.42*vector_c/vector_lscale)); vector_m2
vector_m3 <- pmax(42,ceiling(3.42*vector_c/vector_lscale)); vector_m3
vector_m4 <- pmax(52,ceiling(3.42*vector_c/vector_lscale)); vector_m4
vector_m5 <- pmax(62,ceiling(3.42*vector_c/vector_lscale)); vector_m5
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
setup %>% filter(lscale == 0.05)
setup %>%
  ggplot(aes(x = m, y = HSHS)) + 
  geom_path(col = "black", linewidth = 0.5) +
  geom_path(aes(y = NNHS5, color = "1"), linetype = 1, linewidth = 0.5) +
  geom_path(aes(y = NNHS10, color = "2"), linetype = 1, linewidth = 0.5) +
  geom_path(aes(y = NNHS15, color = "3"), linetype = 1,  linewidth = 0.5) +
  geom_path(aes(y = NNHS20, color = "4"), linetype = 1, linewidth = 0.5) +
  geom_point(aes(y = NNHS5, col = "1"), size = 2, shape = 20) +
  geom_point(aes(y = NNHS10, col = "2"), size =2, shape = 20) +
  geom_point(aes(y = NNHS15, col = "3"), size =2, shape = 20) +
  geom_point(aes(y = NNHS20, col = "4"), size =2, shape = 20) +
  geom_point(aes(fill = "black"), shape = 3) +
  facet_wrap(~lscale, ncol = 4, labeller = label_bquote("\u2113"[k]~"="~.(lscale))) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 5), linetype = "1"), color = "#FC4E07", linewidth = 0.7) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 10), linetype = "2"), color = "#FC4E07", linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 15), linetype = "3"), color = "#FC4E07", linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 20), linetype = "4"), color = "#FC4E07", linewidth = 0.5) +
  scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash"),
                        labels = c("NNNN(m = 5)","NNNN(m = 10)","NNNN(m = 15)","NNNN(m = 20)")) +
  scale_color_manual(values=c("1"="#6fa8dc", "2" = "#3d85c6", "3" = "#476aec", "4" = "#1138cc"), labels = c("NNHS(m = 5)","NNHS(m = 10)","NNHS(m = 15)","NNHS(m = 20)")) +
  scale_fill_manual(values= "black", labels = "HSHS") +
  theme_bw() +
  xlab(bquote(m[d])) +
  ylab(bquote(Computational~Cost)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75,0.25),
        legend.direction = "vertical", 
        legend.box = "horizontal",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.ticks.y = element_blank())
ggsave(filename = "./Rutilities/Computational_Complexity_Plots.png", height = 7, width = 11)


ell <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
c <- pmax(1.20,4.5*ell); c
md <- pmax(22, ceiling(3.42*c/ell)); md
