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

fname <- list.files(path = "./NN10NN10", pattern = "*\\.RData", full.names = FALSE)
nn10nn10_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Model = ) %>% 
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
