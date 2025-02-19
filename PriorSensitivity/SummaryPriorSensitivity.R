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
    mutate(`Elapsed Time (in hours)` = scores_df$`Elapsed Time`/3600) %>%
    select(Method,everything())
  return(fixed_summary)
})


scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  scores_df <- scores_df %>% 
    mutate(ess_tail = min(fixed_summary$ess_tail)) %>%
    mutate(`Elapsed Time` = scores_df$`Elapsed Time`/3600) %>%
    rename(`Elapsed Time (in hours)` = `Elapsed Time`)
  return(scores_df)
})

scores_df <- do.call(rbind, scores_list)
scores_df
scores_df %>% select(Method, ess_tail)


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
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, linetype = Prior), width = 0.10,
                position=position_dodge(0.5), stat="identity") + 
  geom_hline(aes(yintercept = true), linetype = "dotdash", linewidth = 0.5) +
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed, nrow = 2) +
  theme_bw() +
  xlab("") +
  ylab("Posterior 95% CI") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(0.9,0.25))
ggsave(filename = "PriorSensitivity.png", height = 5, width = 11)

library(tidyverse)
fixed_summary %>% group_by(Method,Prior) %>% summarise(`Elapsed Time` = min(`Elapsed Time (in minutes)`)) %>%  arrange(`Elapsed Time`)


scores_df <- do.call(rbind, scores_list)



