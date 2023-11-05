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
  if(exists("scores_df")){
    fixed_summary <- fixed_summary %>% 
      mutate(Method = scores_df$Method) %>%
      mutate(`Elapsed Time (in minutes)` = scores_df$`Elapsed Time`/60) %>%
      select(Method,everything())
    return(fixed_summary)
  }
})


scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("scores_df")){
    return(scores_df)
  }
})

scores_df <- do.call(rbind, scores_list)
scores_df <- scores_df %>% 
  separate(Method,into=c("x","y"),sep="_L") %>%
  mutate(lscale = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  select(-y) %>%
  separate(x,into=c("x1","y1"),sep="_") %>% 
  rename(Method = x1) %>%
  select(-y1)
ggplot(scores_df, aes(x = lscale, y = ES)) + geom_point(aes(col = Method))



## 
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
  mutate(Method = recode(Method, NNNN = 1, NNHS = 2, HSHS = 3)) %>%
  mutate(Prior = recode(Prior, HN = 1, PC = 2, Exp = 3)) 

fixed_summary <- fixed_summary %>% 
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau")))
fixed_summary <- fixed_summary %>% mutate(Prior = factor(Prior, labels = c("Prior Set 1", "Prior Set 2", "Prior Set 3"))) %>% mutate(Method = factor(Method, labels = c("NNNN","NNHS","HSHS")))

ggplot(fixed_summary, aes(x = Method, group = Prior)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, col = Prior), width = 0.25,
                position=position_dodge(0.5), stat="identity") + 
  geom_hline(aes(yintercept = true), linetype = "dotted", linewidth = 0.5) +
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed) +
  theme_bw() +
  xlab("") +
  ylab("Posterior median (95% CI)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = "PriorSensitivity.png", height = 5, width = 9)

library(tidyverse)
fixed_summary %>% group_by(Method,Prior) %>% summarise(`Elapsed Time` = min(`Elapsed Time (in minutes)`)) %>%  arrange(`Elapsed Time`)


