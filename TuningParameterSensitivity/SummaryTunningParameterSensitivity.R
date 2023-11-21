# module load gcc/11.3.0
# module load r/4.2.1
rm(list=ls())
graphics.off()
library(tidyverse)
library(ggplot2)
library(magrittr)
library(coda)
library(fields)
library(lubridate)
library(scoringRules)

fname <- list.files(path = "FullResults", pattern = "^Full_.*\\.RData", full.names = TRUE)
full_fixed_summary_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(Method = "Full") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  cat("i", i,"\n")
  fixed_summary <- fixed_summary %>% select(Method, everything())
  return(fixed_summary)
})
full_fixed_summary_df <- do.call(rbind, full_fixed_summary_list)

full_scores_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(Method = "Full") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
  scores_df <- scores_df %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  scores_df <- scores_df %>% select(Method, everything())
  return(scores_df)
})
full_scores_df <- do.call(rbind, full_scores_list)

## NNNN
fname <- list.files(path = "NNNNResults", pattern = "^NNNN_.*\\.RData", full.names = TRUE)
nnnn_fixed_summary_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(Method = "NNNN") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  fixed_summary <- fixed_summary %>% select(Method, everything())
  return(fixed_summary)
})
nnnn_fixed_summary_df <- do.call(rbind, nnnn_fixed_summary_list)


nnnn_scores_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(Method = "NNNN") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
  scores_df <- scores_df %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  scores_df <- scores_df %>% select(Method, everything())
  return(scores_df)
})
nnnn_scores_df <- do.call(rbind, nnnn_scores_list)

## HSHS1
fname <- list.files(path = "./c1.5HSHS", pattern = "^HSHS_.*\\.RData", full.names = TRUE)
hshs_fixed_summary_list1 <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(m1 = m1, m2 = m2, c = c[1])
  fixed_summary <- fixed_summary %>% mutate(Method = "HSHS") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  fixed_summary <- fixed_summary %>% select(Method, m1,m2, c, everything())
  return(fixed_summary)
})
hshs_fixed_summary_df1 <- do.call(rbind, hshs_fixed_summary_list1)

hshs_scores_list1 <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(m1 = m1, m2 = m2, c = c[1])
  scores_df <- scores_df %>% mutate(Method = "HSHS") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
  scores_df <- scores_df %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  scores_df <- scores_df %>% select(Method, m1,m2, c, everything())
  return(scores_df)
})
hshs_scores_df1 <- do.call(rbind, hshs_scores_list1)

## HSHS2
fname <- list.files(path = "./c1.5HSHS", pattern = "^HSHS_.*\\.RData", full.names = TRUE)
hshs_fixed_summary_list2 <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(m1 = m1, m2 = m2, c = c[1])
  fixed_summary <- fixed_summary %>% mutate(Method = "HSHS") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i])
  fixed_summary <- fixed_summary %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  fixed_summary <- fixed_summary %>% select(Method, m1,m2, c, everything())
  return(fixed_summary)
})
hshs_fixed_summary_df2 <- do.call(rbind, hshs_fixed_summary_list2)

hshs_scores_list2 <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(m1 = m1, m2 = m2, c = c[1])
  scores_df <- scores_df %>% mutate(Method = "HSHS") 
  scores_df <- scores_df %>%  mutate(Setup = fname[i])
  scores_df <- scores_df %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  scores_df <- scores_df %>% select(Method, m1,m2, c, everything())
  return(scores_df)
})
hshs_scores_df2 <- do.call(rbind, hshs_scores_list2)


## HSHSc1p2ell
fname <- list.files(path = "./c1p2ell", pattern = "^HSHS_.*\\.RData", full.names = TRUE)
hshs_fixed_summary_list3 <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(m1 = m1, m2 = m2, c = c[1])
  fixed_summary <- fixed_summary %>% mutate(Method = "HSHS") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  fixed_summary <- fixed_summary %>% select(Method, m1,m2, c, everything())
  return(fixed_summary)
})
hshs_fixed_summary_df3 <- do.call(rbind, hshs_fixed_summary_list3)

hshs_scores_list3 <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(m1 = m1, m2 = m2, c = c[1])
  scores_df <- scores_df %>% mutate(Method = "HSHS") 
  scores_df <- scores_df %>% mutate(lscale = fixed_summary$true[fixed_summary$variable=="ell1"])
  scores_df <- scores_df %>%  mutate(Setup = fname[i])
  scores_df <- scores_df %>% select(Method, m1,m2, c, everything())
  return(scores_df)
})
hshs_scores_df3 <- do.call(rbind, hshs_scores_list3)

#### Merge the results
full_fixed_summary_df <- full_fixed_summary_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>% select(-x) %>% select(Method, Setup, everything())

nnnn_fixed_summary_df <- nnnn_fixed_summary_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs1_fixed_summary_df <- hshs_fixed_summary_df1 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs2_fixed_summary_df <- hshs_fixed_summary_df2 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs3_fixed_summary_df <- hshs_fixed_summary_df3 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())


full_scores_df <- full_scores_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>% select(-x) %>% select(Method, Setup, everything())

nnnn_scores_df <- nnnn_scores_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs1_scores_df <- hshs_scores_df1 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs2_scores_df <- hshs_scores_df2 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs3_scores_df <- hshs_scores_df3 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

lss <- c("full_fixed_summary_df","nnnn_fixed_summary_df","hshs1_fixed_summary_df","hshs2_fixed_summary_df", "hshs3_fixed_summary_df", "full_scores_df","nnnn_scores_df","hshs1_scores_df", "hshs2_scores_df","hshs3_scores_df")
save(list = lss,  file = "Short_Results.rda")

####
rm(list=ls())
load("./TuningParameterSensitivity/Short_Results.rda")

full_fixed_summary_df <-full_fixed_summary_df %>%
  filter(lscale %in% c(0.1,0.4)) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau")))

nnnn_fixed_summary_df %>% 
  filter(lscale %in% c(0.1,0.4)) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) %>%
  ggplot(aes(x = factor(m1))) + 
  geom_point(aes(y = `50%`), size = 1, shape = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1, linewidth = 0.25) + 
  geom_hline(aes(yintercept = true), linetype = "dotted", linewidth = 0.5) +
  ggh4x::facet_grid2(lscale~Pars, scales = "free_y", independent = "y", labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())



hshs1_fixed_summary_df %>% 
  filter(Setup %in% c(5,10,15,20,25)) %>% 
  ggplot(aes(x = m1, y= mean)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  geom_point() + 
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~variable, scales = "free_y") 


hshs1_fixed_summary_df %>% 
  filter(lscale %in% c(0.1,0.4)) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) %>%
  ggplot(aes(x = factor(m1))) + 
  geom_point(aes(y = `50%`), size = 1, shape = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1, linewidth = 0.25) + 
  geom_hline(aes(yintercept = true), linetype = "dotted", linewidth = 0.5) +
  ggh4x::facet_grid2(lscale~Pars, scales = "free_y", independent = "y", labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

hshs2_fixed_summary_df %>% 
  filter(lscale %in% c(0.1,0.4)) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) %>%
  ggplot(aes(x = factor(m1))) + 
  geom_point(aes(y = `50%`), size = 1, shape = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1, linewidth = 0.25) + 
  geom_hline(aes(yintercept = true), linetype = "dotted", linewidth = 0.5) +
  ggh4x::facet_grid2(lscale~Pars, scales = "free_y", independent = "y", labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

hshs3_fixed_summary_df %>% 
  filter(lscale %in% c(0.1,0.4)) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) %>%
  ggplot(aes(x = factor(m1))) + 
  geom_point(aes(y = `50%`), size = 1, shape = 1) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1, linewidth = 0.25) + 
  geom_hline(aes(yintercept = true), linetype = "dotted", linewidth = 0.5) +
  ggh4x::facet_grid2(lscale~Pars, scales = "free_y", independent = "y", labeller = label_parsed) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())
