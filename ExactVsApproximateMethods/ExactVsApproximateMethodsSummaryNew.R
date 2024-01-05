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
fname <- list.files(path = "./Exact", pattern = "*\\.RData", full.names = TRUE)
exact_fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Value = paste0(round(mean,2),"(",round(sd,2),")")) %>%
    select(variable, Value) %>%
    spread(variable, Value) %>%
    select(`theta[1]`, `theta[2]`, `theta[3]`, gamma, sigma1, sigma2, ell1, ell2, tau) %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 0) %>%
    mutate(m = NA, c = NA, m_d = NA) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
exact_fixed_summary <- do.call(rbind, exact_fixed_summary_list)


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


fname <- list.files(path = "./NN15NN15", pattern = "*\\.RData", full.names = TRUE)
nn15nn15_fixed_summary_list <- lapply(1:length(fname), function(node){
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
    mutate(m = 15, c = NA, m_d = NA) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn15nn15_fixed_summary <- do.call(rbind, nn15nn15_fixed_summary_list)


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


fname <- list.files(path = "./NN15HS22", pattern = "*\\.RData", full.names = TRUE)
nn15hs22_fixed_summary_list <- lapply(1:length(fname), function(node){
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
    mutate(m = 15, c = NA, m_d = 22) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn15hs22_fixed_summary <- do.call(rbind, nn15hs22_fixed_summary_list)


fname <- list.files(path = "./NN15HS32", pattern = "*\\.RData", full.names = TRUE)
nn15hs32_fixed_summary_list <- lapply(1:length(fname), function(node){
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
    mutate(m = 15, c = NA, m_d = 32) %>%
    select(DataSet, Method, m, c, m_d, everything())
  return(fixed_summary)
})
nn15hs32_fixed_summary <- do.call(rbind, nn15hs32_fixed_summary_list)


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


fixed_summary <- rbind(exact_fixed_summary, nn10nn10_fixed_summary, nn15nn15_fixed_summary, nn10hs22_fixed_summary, nn10hs32_fixed_summary, nn15hs22_fixed_summary, nn15hs32_fixed_summary, hs22hs22_fixed_summary, hs32hs32_fixed_summary)
fixed_summary %>% arrange(DataSet,Method)

##################################################################################
#### Scoring rules
##################################################################################
#fname <- list.files(path = "./Exact", pattern = "*\\.RData", full.names = TRUE)
#exact_scores_list <- lapply(1:length(fname), function(node){
#  load(fname[node])
#  scores_df <- scores_df %>% 
#    mutate(`ET` = `Elapsed Time`/3600) %>%
#    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
#    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
#    select(-Y) %>%
#    mutate(Method = 0) %>%
#    mutate(m = NA, c = NA, m_d = NA) %>%
#    select(-`Elapsed Time`) %>%
#    select(DataSet, Method, m, c, m_d, everything()) %>%
#    mutate(mESS = round(min(fixed_summary$ess_tail)))
#  return(scores_df)
#})
#exact_scores_df <- do.call(rbind, exact_scores_list)


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


fname <- list.files(path = "./NN15NN15", pattern = "*\\.RData", full.names = TRUE)
nn15nn15_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 15, c = NA, m_d = NA) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn15nn15_scores_df <- do.call(rbind, nn15nn15_scores_list)

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


fname <- list.files(path = "./NN15HS22", pattern = "*\\.RData", full.names = TRUE)
nn15hs22_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 15, c = NA, m_d = 22) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn15hs22_scores_df <- do.call(rbind, nn15hs22_scores_list)


fname <- list.files(path = "./NN15HS32", pattern = "*\\.RData", full.names = TRUE)
nn15hs32_scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(`ET` = `Elapsed Time`/3600) %>%
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 15, c = NA, m_d = 32) %>%
    select(-`Elapsed Time`) %>%
    select(DataSet, Method, m, c, m_d, everything()) %>%
    mutate(mESS = round(min(fixed_summary$ess_tail)))
  return(scores_df)
})
nn15hs32_scores_df <- do.call(rbind, nn15hs32_scores_list)


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

scores_df <- rbind(nn10nn10_scores_df, nn15nn15_scores_df, nn10hs22_scores_df, nn10hs32_scores_df, nn15hs22_scores_df, nn15hs32_scores_df, hs22hs22_scores_df, hs32hs32_scores_df)
scores_df %>% arrange(DataSet,Method) %>% xtable::xtable()

fixed_summary %>% arrange(DataSet,Method) %>% xtable::xtable()



##################################################################################
## KDE
##################################################################################
fname <- list.files(path = "./Exact", pattern = "*\\.RData", full.names = TRUE)
exact_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 0) %>%
    mutate(m = 0)
  return(kde_df)
})
exact_kde_df <- do.call(rbind, exact_kde_list)


fname <- list.files(path = "./NN10NN10", pattern = "*\\.RData", full.names = TRUE)
nn10nn10_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 1)
  return(kde_df)
})
nn10nn10_kde_df <- do.call(rbind, nn10nn10_kde_list)


fname <- list.files(path = "./NN15NN15", pattern = "*\\.RData", full.names = TRUE)
nn15nn15_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 2)
  return(kde_df)
})
nn15nn15_kde_df <- do.call(rbind, nn15nn15_kde_list)

fname <- list.files(path = "./NN10HS22", pattern = "*\\.RData", full.names = TRUE)
nn10hs22_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 3)
  return(kde_df)
})
nn10hs22_kde_df <- do.call(rbind, nn10hs22_kde_list)

fname <- list.files(path = "./NN10HS32", pattern = "*\\.RData", full.names = TRUE)
nn10hs32_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 4)
  return(kde_df)
})
nn10hs32_kde_df <- do.call(rbind, nn10hs32_kde_list)


fname <- list.files(path = "./NN15HS22", pattern = "*\\.RData", full.names = TRUE)
nn15hs22_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 5)
  return(kde_df)
})
nn15hs22_kde_df <- do.call(rbind, nn15hs22_kde_list)


fname <- list.files(path = "./NN15HS32", pattern = "*\\.RData", full.names = TRUE)
nn15hs32_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 6)
  return(kde_df)
})
nn15hs32_kde_df <- do.call(rbind, nn15hs32_kde_list)


fname <- list.files(path = "./HS22HS22", pattern = "*\\.RData", full.names = TRUE)
hs22hs22_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = 7)
  return(kde_df)
})
hs22hs22_kde_df <- do.call(rbind, hs22hs22_kde_list)



fname <- list.files(path = "./HS32HS32", pattern = "*\\.RData", full.names = TRUE)
hs32hs32_kde_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  kde_df <- kde_df %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = 8)
  return(kde_df)
})
hs32hs32_kde_df <- do.call(rbind, hs32hs32_kde_list)


kde_df <- rbind(exact_kde_df, nn10nn10_kde_df, nn15nn15_kde_df, nn10hs22_kde_df, nn10hs32_kde_df, nn15hs22_kde_df, nn15hs32_kde_df, hs22hs22_kde_df, hs32hs32_kde_df)
kde_df <- kde_df %>% mutate(MethodFactor = factor(m , label = c("Exact", "NNNN(10)", "NNNN(15)", "NNHS(10, 22)", "NNHS(10, 32)", "NNHS(15, 22)", "NNHS(15, 32)", "HSHS(22)", "HSHS(32)")))
                         

kde_df %>% 
  filter(DataSet == 1) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = lci, ymax = uci, fill = Key), alpha = 0.35) +
  geom_line(aes(x = x, y = d, col = Key), linetype = "dashed", linewidth = 0.35) +
  facet_wrap(~MethodFactor, ncol = 3) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")
ggsave(filename = "ExactVsApproximateMethodsSummary_gamma_0.75_0.20.png", height = 6, width = 11)


kde_df %>% 
  filter(DataSet == 2) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = lci, ymax = uci, fill = Key), alpha = 0.35) +
  geom_line(aes(x = x, y = d, col = Key), linetype = "dashed", linewidth = 0.35) +
  facet_wrap(~MethodFactor, ncol = 3) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")
ggsave(filename = "ExactVsApproximateMethodsSummary_gamma_0.75_0.50.png", height = 6, width = 11)


kde_df %>% 
  filter(DataSet == 3) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = lci, ymax = uci, fill = Key), alpha = 0.35) +
  geom_line(aes(x = x, y = d, col = Key), linetype = "dashed", linewidth = 0.35) +
  facet_wrap(~MethodFactor, ncol = 3) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")
ggsave(filename = "ExactVsApproximateMethodsSummary_gamma_1.5_0.20.png", height = 6, width = 11)


kde_df %>% 
  filter(DataSet == 4) %>%
  ggplot() +
  geom_ribbon(aes(x = x, ymin = lci, ymax = uci, fill = Key), alpha = 0.35) +
  geom_line(aes(x = x, y = d, col = Key), linetype = "dashed", linewidth = 0.35) +
  facet_wrap(~MethodFactor, ncol = 3) +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")
ggsave(filename = "ExactVsApproximateMethodsSummary_gamma_1.5_0.50.png", height = 6, width = 11)

######################################################################################
# z summary
######################################################################################
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
fname <- list.files(path = "./Exact", pattern = "*\\.RData", full.names = TRUE)
exact_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 0) %>%
    mutate(m = 0) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
exact_z_summary <- do.call(rbind, exact_z_summary_list)
rm(fname)

fname <- list.files(path = "./NN10NN10", pattern = "*\\.RData", full.names = TRUE)
nn10nn10_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 1) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn10nn10_z_summary <- do.call(rbind, nn10nn10_z_summary_list)
rm(fname)

fname <- list.files(path = "./NN15NN15", pattern = "*\\.RData", full.names = TRUE)
nn15nn15_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 1) %>%
    mutate(m = 2) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn15nn15_z_summary <- do.call(rbind, nn15nn15_z_summary_list)


fname <- list.files(path = "./NN10HS22", pattern = "*\\.RData", full.names = TRUE)
nn10hs22_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 3) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn10hs22_z_summary <- do.call(rbind, nn10hs22_z_summary_list)


fname <- list.files(path = "./NN10HS32", pattern = "*\\.RData", full.names = TRUE)
nn10hs32_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 4) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn10hs32_z_summary <- do.call(rbind, nn10hs32_z_summary_list)


fname <- list.files(path = "./NN15HS22", pattern = "*\\.RData", full.names = TRUE)
nn15hs22_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 5) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn15hs22_z_summary <- do.call(rbind, nn15hs22_z_summary_list)


fname <- list.files(path = "./NN15HS32", pattern = "*\\.RData", full.names = TRUE)
nn15hs32_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 2) %>%
    mutate(m = 6) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
nn15hs32_z_summary <- do.call(rbind, nn15hs32_z_summary_list)


fname <- list.files(path = "./HS22HS22", pattern = "*\\.RData", full.names = TRUE)
hs22hs22_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = 7) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
hs22hs22_z_summary <- do.call(rbind, hs22hs22_z_summary_list)

fname <- list.files(path = "./HS32HS32", pattern = "*\\.RData", full.names = TRUE)
hs32hs32_z_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  z_summary <- z_summary %>% 
    mutate(Y = strsplit(fname[node],"_")[[1]][2]) %>%
    mutate(DataSet = as.numeric(gsub(".*?([0-9]+).*", "\\1", Y))) %>%
    select(-Y) %>%
    mutate(Method = 3) %>%
    mutate(m = 8) %>%
    select(DataSet, Method, m, everything()) %>%
    mutate(ID = 1:nrow(.))
  return(z_summary)
})
hs32hs32_z_summary <- do.call(rbind, hs32hs32_z_summary_list)


z_summary <- rbind(
  exact_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn10nn10_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn15nn15_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn10hs22_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn10hs32_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn15hs22_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  nn15hs32_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  hs22hs22_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"),
  hs32hs32_z_summary %>% select("DataSet", "Method", "m", "ID", "post.q2.5", "post.q50", "post.q97.5", "z"))

z_summary %>% arrange(DataSet,Method)

z_summary <- z_summary %>% mutate(MethodFactor = factor(m , label = c("Exact", "NNNN(10)", "NNNN(15)", "NNHS(10, 22)", "NNHS(10, 32)", "NNHS(15, 22)", "NNHS(15, 32)", "HSHS(22)", "HSHS(32)")))


z_summary %>% 
  filter(DataSet == 1) %>%
  filter(ID %in% seq(1,50,1)) %>%
  ggplot(aes(x = ID)) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), alpha = 0.35) +
  geom_point(aes(y = z)) +
  facet_wrap(~MethodFactor, ncol = 3, scales = "free_y") +
  xlab("Latent spatial effect") +
  ylab("Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = "none")
