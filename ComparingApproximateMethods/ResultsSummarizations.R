rm(list=ls())
graphics.off()

library(tidyverse)

dt <- read.csv("./ComparingApproximateMethods/ConvergenceResults.csv")
dt %>% filter(Model == 2) %>% arrange(node)
dt %>% group_by(Model) %>% summarise(Converged = mean(1-converged))
dt %>% group_by(Model) %>% summarise(NotConverged = sum(1-converged))
dt %>% group_by(Model) %>% summarise(Count = length(1-converged))
dt %>% filter(node == 2)
