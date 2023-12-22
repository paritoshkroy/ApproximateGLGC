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

dt <- read.csv("./ComparingApproximateMethods/ConvergenceResults.csv")
str(dt)
dt <- dt %>% distinct()
str(dt)
dt %>% filter(Model == 2) %>% arrange(node)
dt %>% group_by(Model) %>% summarise(`Percent Divergence` = 100*mean(1-converged))
dt %>% group_by(Model) %>% summarise(`Number Divergence` = sum(1-converged))
dt %>% group_by(Model) %>% summarise(`Number Analysis` = length(1-converged))
dt %>% filter(node == 2)

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



