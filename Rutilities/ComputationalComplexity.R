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

set.seed(1492)
ggplot(data.frame(x = seq(500,5000,l=10)), aes(x = x)) +
  geom_function(fun = NNNNComputationalComplexity, args = list(m = 5), col = 1) +
  geom_function(fun = NNNNComputationalComplexity, args = list(m = 10), col = 2) +
  geom_function(fun = NNHSComputationalComplexity, args = list(m = 5, m1 = 22, m2 = 22), col = 3) +
  geom_function(fun = HSHSComputationalComplexity, args = list(m1 = 22, m2 = 22), col = 4) 



ggplot(data.frame(x = seq(5,20,5)), aes(x = x)) +
  geom_function(fun = function(m,n) { 2*n*m^3}, 
                args = list(n = 1000), col = 1) +
  geom_function(fun = function(m,n) {n*m^3 + (n+1)*(m*m)}, 
                args = list(n = 1000), col = 3) +
  geom_function(fun = function(m,n) {2*(n+1)*(m*m) + n}, 
                args = list(n = 1000), col = 4)


NNNNComputationalComplexity(n = 50, m = 15) / HSHSComputationalComplexity(n = 50, m1 = 42, m2 = 42)

NNNNComputationalComplexity(n = 500, m = 5) / HSHSComputationalComplexity(n = 500, m1 = 53, m2 = 23)


vector_lscale <- seq(0.05,0.7,l=14); vector_lscale
vector_c <- round(pmax(4.5*vector_lscale,1.5),2); vector_c
vector_m <- pmax(round(3.42*vector_c/vector_lscale,0),26); vector_m
library(tidyverse)
setup <- data.frame(lscale = vector_lscale, c = vector_c, m = vector_m)
xtable::xtable(setup)

setup <- setup %>%
  mutate(`m = 5_300` = NNNNComputationalComplexity(n = 300, m = 5)) %>%
  mutate(`m = 10_300` = NNNNComputationalComplexity(n = 300, m = 10)) %>%
  mutate(`m = 15_300` = NNNNComputationalComplexity(n = 300, m = 15)) %>%
  mutate(`m = 20_300` = NNNNComputationalComplexity(n = 300, m = 20)) %>%
  mutate(`m = 25_300` = NNNNComputationalComplexity(n = 300, m = 25)) %>%
  mutate(`HSHS_300` = HSHSComputationalComplexity(n = 300, m1 = m, m2 = m))

setup <- setup %>%
  mutate(`m = 5_500` = NNNNComputationalComplexity(n = 500, m = 5)) %>%
  mutate(`m = 10_500` = NNNNComputationalComplexity(n = 500, m = 10)) %>%
  mutate(`m = 15_500` = NNNNComputationalComplexity(n = 500, m = 15)) %>%
  mutate(`m = 20_500` = NNNNComputationalComplexity(n = 500, m = 20)) %>%
  mutate(`m = 25_500` = NNNNComputationalComplexity(n = 500, m = 25)) %>%
  mutate(`HSHS_500` = HSHSComputationalComplexity(n = 500, m1 = m, m2 = m))

setup <- setup %>%
  mutate(`m = 5_1000` = NNNNComputationalComplexity(n = 1000, m = 5)) %>%
  mutate(`m = 10_1000` = NNNNComputationalComplexity(n = 1000, m = 10)) %>%
  mutate(`m = 15_1000` = NNNNComputationalComplexity(n = 1000, m = 15)) %>%
  mutate(`m = 20_1000` = NNNNComputationalComplexity(n = 1000, m = 20)) %>%
  mutate(`m = 25_1000` = NNNNComputationalComplexity(n = 1000, m = 25)) %>%
  mutate(`HSHS_1000` = HSHSComputationalComplexity(n = 1000, m1 = m, m2 = m))

setup <- setup %>%
  mutate(`m = 5_2000` = NNNNComputationalComplexity(n = 2000, m = 5)) %>%
  mutate(`m = 10_2000` = NNNNComputationalComplexity(n = 2000, m = 10)) %>%
  mutate(`m = 15_2000` = NNNNComputationalComplexity(n = 2000, m = 15)) %>%
  mutate(`m = 20_2000` = NNNNComputationalComplexity(n = 2000, m = 20)) %>%
  mutate(`m = 25_2000` = NNNNComputationalComplexity(n = 2000, m = 25)) %>%
  mutate(`HSHS_2000` = HSHSComputationalComplexity(n = 2000, m1 = m, m2 = m))

setup <- setup %>% gather(Key,Value,-lscale,-c,-m)
setup <- setup %>% separate(Key, into = c("Method","n"), sep = "_")
table(setup$Method)

setup %>%
  mutate(n = recode(n, `300`=0, `500`=1, `1000`=2, `2000` =3)) %>%
  mutate(n = factor(n, labels = c("n = 300","n = 500","n = 1000","n = 2000"))) %>%
  mutate(Method = recode(Method, `m = 5`=1, `m = 10`=2, `m = 15` =3, `m = 20` = 4, `m = 25` = 5, `HSHS` = 6)) %>%
  mutate(Method = factor(Method, labels = c("m = 5","m = 10","m = 15", "m = 20", "m = 25", "HSHS"))) %>%
  ggplot(aes(x = lscale, y = Value)) + 
  geom_path(aes(col = Method), linewidth = 0.5, alpha = 2) +
  scale_colour_brewer(palette = "Spectral") +
  facet_wrap(~n, ncol = 2, scales = "free_y") +
  theme_bw() +
  xlab(bquote(Lengthscale)) +
  ylab(bquote(Computational~Complexity)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.ticks.y = element_blank())
ggsave(filename = "Computational_Complexity_Plots.png", height = 6, width = 10)


##
vector_lscale <- seq(0.1,0.5,l=9); vector_lscale
vector_c <- round(pmax(4.5*vector_lscale,1.2),2); vector_c
round(3.42*vector_c/vector_lscale,0)
vector_m1 <- pmax(round(3.42*vector_c/vector_lscale,0),22); vector_m1
vector_m2 <- pmax(round(3.42*vector_c/vector_lscale,0),32); vector_m2
vector_m3 <- pmax(round(3.42*vector_c/vector_lscale,0),44); vector_m3
vector_m4 <- pmax(round(3.42*vector_c/vector_lscale,0),51); vector_m4
vector_m5 <- pmax(round(3.42*vector_c/vector_lscale,0),58); vector_m5
library(tidyverse)
setup <- data.frame(lscale = vector_lscale, c = vector_c, m1 = vector_m1, m2 = vector_m2, m3 = vector_m3, m4 = vector_m4, m5 = vector_m5)

xtable::xtable(setup)

setup <- setup %>%
  mutate(`HSHS(22)` = HSHSComputationalComplexity(n = 1000, m1 = m1, m2 = m1)) %>%
  mutate(`HSHS(32)` = HSHSComputationalComplexity(n = 1000, m1 = m2, m2 = m2)) %>%
  mutate(`HSHS(44)` = HSHSComputationalComplexity(n = 1000, m1 = m3, m2 = m3)) %>%
  mutate(`HSHS(51)` = HSHSComputationalComplexity(n = 1000, m1 = m4, m2 = m4)) %>%
  mutate(`HSHS(58)` = HSHSComputationalComplexity(n = 1000, m1 = m5, m2 = m5))

setup
setup <- setup %>% select(-c,-m1,-m2,-m3,-m4,-m5) %>% gather(Method,Value,-lscale)

table(setup$Method)

setup %>%
  mutate(Method = recode(Method, `HSHS(22)` = 1, `HSHS(32)` = 2, `HSHS(44)` = 3, `HSHS(51)` = 4, `HSHS(58)` = 5)) %>%
  mutate(Method = factor(Method, labels = c("HSHS(22)", "HSHS(32)","HSHS(44)","HSHS(51)","HSHS(58)"))) %>%
  ggplot(aes(x = lscale, y = Value)) + 
  geom_path(linewidth = 0.5, linetype = "solid") +
  geom_point(size = 1, shape = 1) +
  facet_wrap(~Method, ncol = 3) +
  scale_colour_brewer(palette = "Spectral") +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 5), linetype = "1")) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 10), linetype = "2")) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 15), linetype = "3")) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 20), linetype = "4")) +
  scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash"),
                        labels = c("NNNN(5)","NNNN(10)","NNNN(15)","NNNN(20)")) +
  theme_bw() +
  xlab(bquote(Lengthscale)) +
  ylab(bquote(Computational~Complexity)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85,0.25),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.ticks.y = element_blank())
ggsave(filename = "Computational_Complexity_Plots.png", height = 5, width = 10)


###
## c = 1.5
vector_lscale <- seq(0.05,0.5,l=10); vector_lscale
vector_c <- round(pmax(4.5*vector_lscale,1.2),2); vector_c
round(3.42*vector_c/vector_lscale,0)
vector_m1 <- pmax(22,round(3.42*vector_c/vector_lscale,0)); vector_m1
vector_m2 <- pmax(28,round(3.42*vector_c/vector_lscale,0)); vector_m2
vector_m3 <- pmax(35,round(3.42*vector_c/vector_lscale,0)); vector_m3
vector_m4 <- pmax(42,round(3.42*vector_c/vector_lscale,0)); vector_m4
vector_m5 <- pmax(52,round(3.42*vector_c/vector_lscale,0)); vector_m5
vector_m6 <- pmax(60,round(3.42*vector_c/vector_lscale,0)); vector_m6
vector_m7 <- pmax(82,round(3.42*vector_c/vector_lscale,0)); vector_m7

#vector_m1 <- 22; vector_m1
#vector_m2 <- 26; vector_m2
#vector_m3 <- 32; vector_m3
#vector_m4 <- 44; vector_m4
#vector_m5 <- 51; vector_m5
#vector_m6 <- 58; vector_m6
library(tidyverse)
setup <- tibble(lscale = vector_lscale, c = vector_c, m1 = vector_m1, m2 = vector_m2, m3 = vector_m3, m4 = vector_m4, m5 = vector_m5, m6 = vector_m6, m7 = vector_m7)
setup <- setup %>% gather(x,m,-c,-lscale) %>% select(-x) %>% distinct()
setup
table(setup$m)

setup <- setup %>% mutate(cs = HSHSComputationalComplexity(n = 1000, m1 = m, m2 = m)) 
setup
setup <- setup %>% mutate(Method = "HSHS")
setup
setup %>% filter(lscale == 0.10)
setup %>%
  ggplot(aes(x = factor(m), y = cs)) + 
  geom_path(aes(group = 1), linewidth = 0.7) +
  geom_point(size = 2, shape = 1) +
  facet_wrap(~lscale, ncol = 4, labeller = label_bquote("\u2113"[k]~"="~.(lscale))) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 5), linetype = "1"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 10), linetype = "2"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 15), linetype = "3"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 20), linetype = "4"), linewidth = 0.5) +
  scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash"),
                        labels = c("NNNN(5)","NNNN(10)","NNNN(15)","NNNN(20)")) +
  theme_bw() +
  xlab(bquote(Number~of~basis~functions~under~HSHS)) +
  ylab(bquote(Computational~Cost)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75,0.125),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.ticks.y = element_blank())
ggsave(filename = "Computational_Complexity_Plots.png", height = 7, width = 11)


## c = 1.30
vector_lscale <- seq(0.05,0.5,l=10); vector_lscale
round(4.5*vector_lscale,2)
vector_c <- pmax(1.30, round(4.5*vector_lscale,2)); vector_c
vector_c <- 1.2 + 2*vector_lscale; vector_c
cbind(vector_lscale, vector_c, round(3.42*vector_c/vector_lscale,0))
vector_m1 <- pmax(23,round(3.42*vector_c/vector_lscale,0)); vector_m1
vector_m2 <- pmax(27,round(3.42*vector_c/vector_lscale,0)); vector_m2
vector_m3 <- pmax(30,round(3.42*vector_c/vector_lscale,0)); vector_m3
vector_m4 <- pmax(34,round(3.42*vector_c/vector_lscale,0)); vector_m4
vector_m5 <- pmax(48,round(3.42*vector_c/vector_lscale,0)); vector_m5
vector_m6 <- pmax(89,round(3.42*vector_c/vector_lscale,0)); vector_m6

library(tidyverse)
setup <- tibble(lscale = vector_lscale, c = vector_c, m1 = vector_m1, m2 = vector_m2, m3 = vector_m3, m4 = vector_m4, m5 = vector_m5, m6 = vector_m6)
setup <- setup %>% gather(x,m,-c,-lscale) %>% select(-x) %>% distinct()
setup
table(setup$m)

setup <- setup %>% mutate(cs = HSHSComputationalComplexity(n = 1000, m1 = m, m2 = m)) 
setup
setup <- setup %>% mutate(Method = "HSHS")
setup
setup %>% filter(lscale == 0.10)
setup %>%
  ggplot(aes(x = m, y = cs)) + 
  geom_path(linewidth = 0.7) +
  geom_point(size = 2, shape = 1) +
  facet_wrap(~lscale, ncol = 4, labeller = label_bquote("\u2113"[k]~"="~.(lscale))) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 5), linetype = "1"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 10), linetype = "2"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 15), linetype = "3"), linewidth = 0.5) +
  geom_hline(aes(yintercept = NNNNComputationalComplexity(n = 1000, m = 20), linetype = "4"), linewidth = 0.5) +
  scale_linetype_manual(values=c("dotted","dotdash","dashed","longdash"),
                        labels = c("NNNN(5)","NNNN(10)","NNNN(15)","NNNN(20)")) +
  theme_bw() +
  xlab(bquote(Number~of~basis~functions~under~HSHS)) +
  ylab(bquote(Computational~Cost)) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.75,0.125),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.ticks.y = element_blank())
ggsave(filename = "Computational_Complexity_Plots.png", height = 7, width = 11)
