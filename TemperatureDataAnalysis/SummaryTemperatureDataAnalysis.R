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
fname <- list.files(pattern = "*\\.RData", full.names = FALSE)

fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Model = fname[node]) %>% 
    select(Model,everything())
  return(fixed_summary)
})
fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary <- fixed_summary %>% separate(Model, into = c("x","y"), sep = "_G") %>%
  mutate(x = recode(x, Full="Full0_Full0")) %>% select(-y) %>% rename(Model=x) %>%
  separate(Model, into = c("x","y"), sep = "_") %>%
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  mutate(Method = recode(Method, FullFull="Full")) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,everything())

fixed_summary <- fixed_summary %>% mutate(MethodCharc = paste0(Method,"(",m1,",",m2,")")) 
fixed_summary <- fixed_summary %>% 
  mutate(MethodFactor = recode(MethodCharc, "NNNN(6,6)" = 1, "NNNN(10,10)" = 2, "NNNN(15,15)" = 3, "NNNN(20,20)" = 4, "HSHS(22,22)" = 5, "HSHS(30,30)" = 6, "HSHS(36,36)" = 7, "HSHS(42,42)" = 8))
fixed_summary <- fixed_summary %>% 
  mutate(MethodFactor = factor(MethodFactor, labels = c("NNNN(6,6)","NNNN(10,10)","NNNN(15,15)", "NNNN(20,20)", "HSHS(22,22)", "HSHS(30,30)", "HSHS(36,36)", "HSHS(42,42)")))

fixed_summary <- fixed_summary %>% 
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) 

ggplot(fixed_summary, aes(x = MethodFactor)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25,
                position=position_dodge(0.5), stat="identity") + 
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed, nrow = 3) +
  theme_bw() +
  xlab("") +
  ylab("Posterior median (95% CI)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = "FixedParsSummaries.png", height = 5, width = 9)

###################################################################
## Predictive Scores
###################################################################
scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("scores_df")) return(scores_df %>% mutate(Model = fname[node]))
})
scores_df <- do.call(rbind, scores_list)
scores_df <- scores_df %>% separate(Model, into = c("x","y"), sep = "_G") %>%
  mutate(x = recode(x, Full="Full0_Full0")) %>% select(-y) %>% rename(Model=x) %>%
  separate(Model, into = c("x","y"), sep = "_") %>%
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  mutate(Method = recode(Method, FullFull="Full")) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,everything())

scores_df <- scores_df %>% mutate(MethodCharc = paste0(Method,"(",m1,",",m2,")")) 
scores_df <- scores_df %>% 
  mutate(MethodFactor = recode(MethodCharc, "NNNN(6,6)" = 1, "NNNN(10,10)" = 2, "NNNN(15,15)" = 3, "NNNN(20,20)" = 4, "HSHS(22,22)" = 5, "HSHS(30,30)" = 6, "HSHS(36,36)" = 7, "HSHS(42,42)" = 8))
scores_df <- scores_df %>% 
  mutate(MethodFactor = factor(MethodFactor, labels = c("NNNN(6,6)","NNNN(10,10)","NNNN(15,15)", "NNNN(20,20)", "HSHS(22,22)", "HSHS(30,30)", "HSHS(36,36)", "HSHS(42,42)")))
scores_df %>% select(MethodFactor, MAE, RMSE, CVG, CRPS, IS, ES, logs, `Elapsed Time`) %>% arrange(MethodFactor) %>% xtable::xtable()


