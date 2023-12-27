# module load gcc/11.3.0
# module load r/4.2.1
###########################################################################
### To the Compute Canada
###########################################################################
#module load gcc/11.3.0
#module load r/4.2.1
rm(list=ls())
graphics.off()
library(tidyverse)
library(ggplot2)
library(magrittr)
library(coda)
library(fields)
library(lubridate)
library(scoringRules)

fname <- list.files(path = ".", pattern = "*\\.RData", full.names = TRUE)
#fname <- list.files(path = ".", pattern = "^NNNN_.*\\.RData", full.names = TRUE)

fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("fixed_summary")){
  fixed_summary <- fixed_summary %>% 
    mutate(Method = scores_df$Method) %>%
    mutate(`Elapsed Time (in hours)` = scores_df$`Elapsed Time`/3600) %>%
    mutate(node =  as.numeric(gsub(".*?([0-9]+).*", "\\1", fname[node]))) %>%
    select(Method,everything())
  return(fixed_summary)
  }
})

fixed_summary_df <- do.call(rbind, fixed_summary_list)

scores_df_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("scores_df")){
  scores_df <- scores_df %>% 
    mutate(`Elapsed Time (in hours)` = scores_df$`Elapsed Time`/3600) %>%
    mutate(node =  as.numeric(gsub(".*?([0-9]+).*", "\\1", fname[node]))) %>%
    select(Method,everything())
  return(scores_df)
  }
})
scores_df <- do.call(rbind, scores_df_list)

save(fixed_summary_df, scores_df, file = "ResultsSummarizations.RData")

############################################################################
# Local PC
############################################################################
rm(list=ls())
graphics.off()
library(tidyverse)
library(ggplot2)
library(magrittr)
library(coda)
library(fields)
library(lubridate)
library(scoringRules)

dt <- read.csv("./ComparingApproximateMethods/ConvergenceResults.csv")
str(dt)
dt <- dt %>% distinct()
str(dt)
dt %>% filter(Model == 2) %>% arrange(node)
dt %>% group_by(Model) %>% summarise(`Percent Divergence` = 100*mean(1-converged))
dt %>% group_by(Model) %>% summarise(`Number Divergence` = sum(1-converged))
dt %>% group_by(Model) %>% summarise(`Number Analysis` = length(1-converged))
dt %>% filter(node == 2)

head(dt)
dt_converged <- dt %>% filter(converged==1)
head(dt_converged)
Model1 <- dt_converged %>% filter(Model == 1) %>% rename(Model1=Model)
Model2 <- dt_converged %>% filter(Model == 2) %>% rename(Model2=Model)
Model3 <- dt_converged %>% filter(Model == 3) %>% rename(Model3=Model)

dt_all_converged <- inner_join(inner_join(Model1, Model2, by = c("node","converged")), Model3, by = c("node","converged"))
head(dt_all_converged)
dt_all_converged <- dt_all_converged %>% gather(Mode,Model,-node,-converged)
dt_all_converged <- dt_all_converged %>% select(-Mode)
dt_all_converged <- dt_all_converged %>% mutate(Model = factor(Model, labels = c("NNNN", "NNHS", "HSHS")))

head(dt_all_converged)
dt_all_converged <- dt_all_converged %>% arrange(node,Model) %>% group_by(node) %>% mutate(nodeID = cur_group_id())

load("./ComparingApproximateMethods/ResultsSummarizations.RData")

head(fixed_summary_df)
table(fixed_summary_df$Method)
fixed_summary_df <- fixed_summary_df %>% 
  mutate(Model = recode(Method, "NNNN_GLGC" = 1, "NNHS_GLGC" = 2, "HSHS_GLGC" = 3)) %>%
  mutate(Model = factor(Model, labels = c("NNNN", "NNHS", "HSHS")))
table(fixed_summary_df$Model)
fixed_summary_df <- inner_join(fixed_summary_df, dt_all_converged, by = c("Model","node"))

fixed_summary_df %>% group_by(Model) %>% summarise(`Elapsed Time (in hours)` = mean(`Elapsed Time (in hours)`))
head(fixed_summary_df)

fixed_summary_gamma <- fixed_summary_df %>% 
  filter(Model == "NNHS" & variable == "gamma") %>% 
  mutate(`2.5%` = `2.5%`/0.4, `97.5%` = `97.5%`/0.4)
fixed_summary_minus_gamma <- fixed_summary_df %>% 
  filter(!(Model == "NNHS" & variable == "gamma"))

rbind(fixed_summary_gamma, fixed_summary_minus_gamma) %>% 
  mutate(Model = factor(Model, labels = c("NNNN", "NNHS", "HSHS"))) %>%
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) %>% 
  ggplot(aes(x = factor(nodeID), group = Model, col = Model)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                position=position_dodge(width=0.5), width = 0.1) +
  geom_hline(aes(yintercept = true), linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed, ncol = 2) +
  ylab("Posterior 95% CI") +
  xlab(bquote(Sample)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        legend.position = c(0.8,0.08),
        legend.direction="horizontal",
        legend.title = element_blank())
ggsave(filename = "./ComparingApproximateMethods/replicated_fixed_parameters.png", height = 9, width = 11)

mESS_dt <- rbind(fixed_summary_gamma, fixed_summary_minus_gamma) %>% 
  mutate(Model = factor(Model, labels = c("NNNN", "NNHS", "HSHS"))) %>%
  group_by(Model, node) %>%
  summarize(mESS = min(ess_tail)) %>%
  ungroup()

names(scores_df)
head(scores_df)

scores_all_df <- scores_df %>% 
  mutate(Model = recode(Method, "NNNN_GLGC" = 1, "NNHS_GLGC" = 2, "HSHS_GLGC" = 3)) %>%
  mutate(Model = factor(Model, labels = c("NNNN", "NNHS", "HSHS"))) %>% 
  select(-"Elapsed Time",-"Method")

scores_all_df <- inner_join(inner_join(scores_all_df, dt_all_converged, by = c("Model","node")), mESS_dt, by = c("Model","node")) %>% select(-node)

scores_df_long <- scores_all_df %>% select(-converged, -node, -nodeID) %>%
  gather(Key, Value, -Model)
  
  
scores_df_long <- scores_df_long %>% 
  mutate(KeyRecode = recode(Key, "MAE" = 1, "RMSE" = 2, "CVG" = 3, "CRPS" = 4, "IS" = 5, "ES" = 6, "logs" = 7, "Elapsed Time (in hours)" = 8, "mESS" = 9)) %>%
  mutate(KeyFactor = factor(KeyRecode, label = c("MAE", "RMSE", "CVG", "CRPS", "IS", "ES", "logs", "Elapsed Time (in hours)", "mESS")))

scores_df_long %>%
  ggplot(aes(x = Model)) +
  geom_boxplot(aes(y = Value)) +
  facet_wrap(~KeyFactor, scales = "free_y") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10))

ggsave(filename = "./ComparingApproximateMethods/replicated_scoring_rules.png", height = 6, width = 11)
