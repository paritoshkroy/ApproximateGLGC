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

### Fixed Parameters
fname <- list.files(path = "./NN10NN10", pattern = "*\\.RData", full.names = TRUE)
nn10nn10_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 10, c = NA, m_d = NA) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn10nn10_fixed_summary <- do.call(rbind, nn10nn10_fixed_summary_list)


fname <- list.files(path = "./NN10HS22", pattern = "*\\.RData", full.names = TRUE)
nn10hs22_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 10, c = NA, m_d = 22) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn10hs22_fixed_summary <- do.call(rbind, nn10hs22_fixed_summary_list)


fname <- list.files(path = "./NN10HS32", pattern = "*\\.RData", full.names = TRUE)
nn10hs32_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 10, c = NA, m_d = 32) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn10hs32_fixed_summary <- do.call(rbind, nn10hs32_fixed_summary_list)

fname <- list.files(path = "./HS22HS22", pattern = "*\\.RData", full.names = TRUE)
hs22hs22_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = NA, c = NA, m_d = 22) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
hs22hs22_fixed_summary <- do.call(rbind, hs22hs22_fixed_summary_list)

fname <- list.files(path = "./HS32HS32", pattern = "*\\.RData", full.names = TRUE)
hs32hs32_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = NA, c = NA, m_d = 32) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
hs32hs32_fixed_summary <- do.call(rbind, hs32hs32_fixed_summary_list)


fixed_summary <- rbind(nn10nn10_fixed_summary, nn10hs22_fixed_summary, nn10hs32_fixed_summary, hs22hs22_fixed_summary, hs32hs32_fixed_summary)
fixed_summary %>% arrange(DataSet,Method)

##################################################################################
#### Scoring rules
##################################################################################

fname <- list.files(path = "./NN10NN10", pattern = "*\\.RData", full.names = TRUE)
nn10nn10_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 10, c = NA, m_d = NA) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn10nn10_scores_df <- do.call(rbind, nn10nn10_scores_list)


fname <- list.files(path = "./NN10HS22", pattern = "*\\.RData", full.names = TRUE)
nn10hs22_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 10, c = NA, m_d = 22) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn10hs22_scores_df <- do.call(rbind, nn10hs22_scores_list)


fname <- list.files(path = "./NN10HS32", pattern = "*\\.RData", full.names = TRUE)
nn10hs32_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 10, c = NA, m_d = 32) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn10hs32_scores_df <- do.call(rbind, nn10hs32_scores_list)


fname <- list.files(path = "./HS22HS22", pattern = "*\\.RData", full.names = TRUE)
hs22hs22_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = NA, c = NA, m_d = 22) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
hs22hs22_scores_df <- do.call(rbind, hs22hs22_scores_list)


fname <- list.files(path = "./HS32HS32", pattern = "*\\.RData", full.names = TRUE)
hs32hs32_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = NA, c = NA, m_d = 32) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
hs32hs32_scores_df <- do.call(rbind, hs32hs32_scores_list)

scores_df <- rbind(nn10nn10_scores_df, nn10hs22_scores_df, nn10hs32_scores_df, hs22hs22_scores_df, hs32hs32_scores_df)
scores_df %>% arrange(DataSet,Method)

fixed_summary %>% arrange(DataSet,Method)
