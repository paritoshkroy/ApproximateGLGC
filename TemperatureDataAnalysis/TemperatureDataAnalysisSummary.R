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

fname <- list.files(path = ".", pattern = "*\\.RData", full.names = FALSE)
fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Model = fname[node]) %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(Model, variable, Value)
  return(fixed_summary)
})
fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary <- fixed_summary %>% 
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9, sigma = 6, ell = 8)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau")))
fixed_summary <- fixed_summary %>% 
  mutate(Method = recode(Model, `NNGP_Temps.RData` = 1, `logNNGP_Temps.RData` = 2, `NNNN_GLGC_Temps.RData` = 3, `NNHS1_GLGC_Temps.RData` = 4, `NNHS2_GLGC_Temps.RData` = 5, `HSHS1_GLGC_Temps.RData` = 6, `HSHS2_GLGC_Temps.RData` = 7)) %>%
  mutate(Method = factor(Method, labels = c("NNGP","log NNGP", "NNNN", "NNHS1", "NNHS2", "HSHS1", "HSHS2"))) %>%
  select(Method,Pars,Value)
fixed_summary %>% spread(Pars,Value) %>% xtable::xtable()

### Scoring rules

fname <- list.files(path = ".", pattern = "*\\.RData", full.names = FALSE)
scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail))) %>%
    select(- `Elapsed Time`)
  return(scores_df)
})
scores_df <- do.call(rbind, scores_list)
scores_df <- scores_df %>% 
  mutate(Method = recode(Method, `NNGP` = 1, `logNNGP` = 2, `NNNN_GLGC` = 3, `NNHS1_GLGC` = 4, `NNHS2_GLGC` = 5, `HSHS1_GLGC` = 6, `HSHS2_GLGC` = 7)) %>%
  mutate(Method = factor(Method, labels = c("NNGP","log NNGP", "NNNN", "NNHS1", "NNHS2", "HSHS1", "HSHS2")))
scores_df %>% arrange(Method) %>% xtable::xtable()


