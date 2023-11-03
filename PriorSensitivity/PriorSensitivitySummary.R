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
###########################################################################
### Gaussian Model
###########################################################################
fname <- list.files(pattern = "*\\.RData", full.names = FALSE)
fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Method = scores_df$Method) %>%
    mutate(`Elapsed Time (in minutes)` = scores_df$`Elapsed Time`/60) %>%
    select(Method,everything())
  return(fixed_summary)
})
fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary %>% filter(variable %in% "ell1")
fixed_summary %>% filter(variable %in% "ell2") 
fixed_summary %>% filter(variable %in% "gamma")
fixed_summary <- fixed_summary %>% 
  separate(Method,into=c("x","y"),sep="_G") %>% 
  rename(Method = x) %>% 
  separate(y, into = c("x","y"), sep = "_") %>% 
  rename(Prior = y) %>% 
  select(-x) %>% 
  mutate(Method = recode(Method, Full = 1, NNNN = 2, NNHS = 3, HSHS = 4)) %>%
  mutate(Prior = recode(Prior, HN = 1, PC = 2, Exp = 3)) 

ggplot(fixed_summary, aes(x = Prior)) + 
  geom_errorbar(aes(ymin=`2.5%`,ymax=`97.5%`)) + 
  geom_point(aes(y = `50%`)) + 
  facet_wrap(Method~variable)

