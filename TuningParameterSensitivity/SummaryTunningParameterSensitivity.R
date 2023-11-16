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

fname <- list.files(path = ".", pattern = "^Full_.*\\.RData", full.names = FALSE)
full_fixed_summary_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(Method = "Full") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  cat("i", i,"\n")
  fixed_summary <- fixed_summary %>% select(Method, everything())
  return(fixed_summary)
})
full_fixed_summary_df <- do.call(rbind, full_fixed_summary_list)

full_scores_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(Method = "Full") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
  scores_df <- scores_df %>% select(Method, everything())
  return(scores_df)
})
full_scores_df <- do.call(rbind, full_scores_list)

## NNNN
fname <- list.files(path = ".", pattern = "^NNNN_.*\\.RData", full.names = FALSE)
nnnn_fixed_summary_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(Method = "NNNN") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% select(Method, everything())
  return(fixed_summary)
})
nnnn_fixed_summary_df <- do.call(rbind, nnnn_fixed_summary_list)


nnnn_scores_list <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(Method = "NNNN") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
  scores_df <- scores_df %>% select(Method, everything())
  return(scores_df)
})
nnnn_scores_df <- do.call(rbind, nnnn_scores_list)

## HSHS1
fname <- list.files(path = ".", pattern = "^HSHS_.*\\.RData", full.names = FALSE)
hshs_fixed_summary_list1 <- lapply(1:length(fname), function(i){
  load(fname[i])
  fixed_summary <- fixed_summary %>% mutate(m1 = m1, m2 = m2, c = c[1])
  fixed_summary <- fixed_summary %>% mutate(Method = "HSHS") 
  fixed_summary <- fixed_summary %>% mutate(Setup = fname[i]) 
  fixed_summary <- fixed_summary %>% select(Method, m1,m2, c, everything())
  return(fixed_summary)
})
hshs_fixed_summary_df1 <- do.call(rbind, hshs_fixed_summary_list1)

hshs_scores_list1 <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(m1 = m1, m2 = m2, c = c[1])
  scores_df <- scores_df %>% mutate(Method = "HSHS") 
  scores_df <- scores_df %>% mutate(Setup = fname[i]) 
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
  fixed_summary <- fixed_summary %>% select(Method, m1,m2, c, everything())
  return(fixed_summary)
})
hshs_fixed_summary_df2 <- do.call(rbind, hshs_fixed_summary_list2)

hshs_scores_list2 <- lapply(1:length(fname), function(i){
  load(fname[i])
  scores_df <- scores_df %>% mutate(m1 = m1, m2 = m2, c = c[1])
  scores_df <- scores_df %>% mutate(Method = "HSHS") 
  scores_df <- scores_df %>%  mutate(Setup = fname[i])
  scores_df <- scores_df %>% select(Method, m1,m2, c, everything())
  return(scores_df)
})
hshs_scores_df2 <- do.call(rbind, hshs_scores_list2)

#### Merge the results

full_fixed_summary_df <- full_fixed_summary_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>% select(-x) %>% select(Method, Setup, everything())

nnnn_fixed_summary_df <- nnnn_fixed_summary_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs1_fixed_summary_df <- hshs_fixed_summary_df1 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs2_fixed_summary_df <- hshs_fixed_summary_df2 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())



full_scores_df <- full_scores_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>% select(-x) %>% select(Method, Setup, everything())

nnnn_scores_df <- nnnn_scores_df %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs1_scores_df <- hshs_scores_df1 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())

hshs2_scores_df <- hshs_scores_df2 %>% separate(Setup, into = c("x","y"), sep = "_") %>% mutate(Setup = as.numeric(gsub(".*?([0-9]+).*", "\\1", y)))  %>% select(-x) %>% select(Method, Setup, everything())


lss <- c("full_fixed_summary_df","nnnn_fixed_summary_df","hshs1_fixed_summary_df","hshs2_fixed_summary_df","full_scores_df","nnnn_scores_df","hshs1_scores_df",
     "hshs2_scores_df")
save(list = lss,  file = "Short_Results.rda")

load("Short_Results.rda")

nnnn_fixed_summary_df %>% 
  filter(Setup %in% c(1,6,11)) %>% 
  ggplot(aes(x = Setup, y= mean)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  geom_point() + 
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~variable, scales = "free_y") 


hshs1_fixed_summary_df %>% 
  filter(Setup %in% c(5,10,15,20,25)) %>% 
  ggplot(aes(x = m1, y= mean)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  geom_point() + 
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~variable, scales = "free_y") 


hshs1_fixed_summary_df %>% 
  filter(Setup %in% c(1,6,11,16,21)) %>% 
  ggplot(aes(x = m1, y= mean)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  geom_point() + 
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~variable, scales = "free_y") 

hshs2_fixed_summary_df %>% 
  filter(Setup %in% c(1,6,11,16,21)) %>% 
  ggplot(aes(x = m1, y= mean)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
  geom_point() + 
  geom_hline(aes(yintercept = true)) +
  facet_wrap(~variable, scales = "free_y") 


ggplot(data = hshs2_fixed_summary_df, aes(x = date)) + 
  geom_path(aes(y = y), col = "dimgray", linewidth = 0.5) + 
  facet_wrap(~varid, ncol = 1, labeller = label_bquote(y[.(varid)])) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))
