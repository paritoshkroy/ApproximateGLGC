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
    mutate(LS = unlist(strsplit(fname[node], "_"))[3]) %>%
    mutate(LS = as.numeric(gsub(".*?([0-9]+).*", "\\1", LS))) %>%
    mutate(Method = paste(unlist(strsplit(fname[node], "_"))[1], unlist(strsplit(fname[node], "_"))[2], sep="_")) %>%
    select(Method,LS,everything())
  return(fixed_summary)
})
fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary %>% filter(variable %in% "ell1") %>% arrange(LS)
fixed_summary %>% filter(variable %in% "ell2") %>% arrange(LS)
fixed_summary %>% filter(variable %in% "gamma") %>% arrange(LS)
fixed_summary %>% arrange(LS == 1)

fixed_summary %>% 
  filter(LS == 1) %>%
  filter(variable %in% "ell1") %>% 
  ggplot(aes(x = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
  geom_point(aes(y = `50%`)) +
  facet_wrap(~LS, nrow = 2, labeller = label_bquote("\u2113"[1]~"="~.(LS))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))

fixed_summary %>% 
  filter(variable %in% "ell2") %>% 
  ggplot(aes(x = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
  geom_point(aes(y = `50%`)) +
  facet_wrap(~LS, nrow = 2, labeller = label_bquote("\u2113"[2]~"="~.(LS))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))

fixed_summary %>% 
  filter(variable %in% "gamma") %>% 
  filter(LS == 2) %>%
  ggplot(aes(x = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
  geom_point(aes(y = `50%`)) +
  facet_wrap(~LS, nrow = 2, labeller = label_bquote(gamma~"="~.(LS))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))


fixed_summary %>% 
  filter(variable %in% "theta[1]") %>% 
  filter(LS == 1) %>%
  ggplot(aes(x = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
  geom_point(aes(y = `50%`)) +
  facet_wrap(~LS, nrow = 2, labeller = label_bquote(gamma~"="~.(LS))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))

#######################
scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("scores_df")) return(scores_df %>% mutate(LS = unlist(strsplit(fname[node], "_"))[[3]]))
})
scores_df <- do.call(rbind, scores_list)
scores_df <- scores_df %>% 
  mutate(LS = as.numeric(gsub(".*?([0-9]+).*", "\\1", LS)))
scores_df <- scores_df %>% separate(Method, into = c("x","y"), sep = "_L")
scores_df <- scores_df %>% rename(Method = x) %>% select(-y)
scores_df <- scores_df %>% separate(Method, into = c("x","y"), sep = "_")
scores_df <- scores_df %>% 
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>% 
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,LS,everything()) %>%
  replace_na(list(m1 = 0, m2 = 0))
scores_df <- scores_df %>% filter(!(m1 == 25))
library(gridExtra)
grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN") & LS >1) %>% 
    ggplot(aes(x = LS, y = ES)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("NNNN") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS") & LS >1) %>% 
    ggplot(aes(x = LS, y = ES)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("NNHS") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","HSHS") & LS >1) %>% 
    ggplot(aes(x = LS, y = ES)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("HSHS") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)), 
  
  ncol = 3)

#### Logs
grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN") & LS >1) %>% 
    ggplot(aes(x = LS, y = logs)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("NNNN") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS") & LS >1) %>% 
    ggplot(aes(x = LS, y = logs)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("NNHS") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","HSHS") & LS >1) %>% 
    ggplot(aes(x = LS, y = logs)) + 
    geom_line(aes(col = factor(m1))) + 
    ggtitle("HSHS") +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.8), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)), 
  
  ncol = 3)

## Following indicates that to estimate ell of minimum 0.2 NN10NN10, NN10HS36 and HS36HS36 provide the same prediction accuracy
scores_df %>% select(-RMSE) %>% mutate(`Elapsed Time` = `Elapsed Time`/3600) %>% filter(LS==2) %>% arrange(`Elapsed Time`)
scores_df %>% select(-RMSE) %>% mutate(`Elapsed Time` = `Elapsed Time`/3600) %>% filter(LS==5) %>% arrange(`Elapsed Time`) %>% print(n=100)
scores_df %>% select(-RMSE) %>% mutate(`Elapsed Time` = `Elapsed Time`/3600) %>% filter(LS==9) %>% arrange(`Elapsed Time`) %>% print(n=100)


fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary %>% group_by(node) %>% summarize(max = max(rhat)) %>% print(n = 33)
fixed_summary %>% filter(variable == "gamma") %>% filter(LS==2)

### CPOP
cpop <- function(x,y){
  if(median(x)<y){
    mean(x>y)
  } else {
    mean(x<y)
  }
}

cpop_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  size_post_samples <- dim(out$post_ypred)[1]
  TT <- dim(out$post_ypred)[2]
  q <- dim(out$post_ypred)[3]
  cpop_summary <- matrix(nrow = TT, ncol = q)
  for(i in 1:TT){
    cpop_summary[i,] <- cbind(cpop(out$post_ypred[,i,1],out$y1_test[i]), cpop(out$post_ypred[,i,2],out$y2_test[i]))
  }
  cpop_summary_df <- as_tibble(matrix2df(cpop_summary, rowindex = rep(unique(out$mypred$date), colindex = 1:2))) %>% rename(date = rindex, varid = cindex, cpop = x)
  out <- inner_join(out$mypred, cpop_summary_df, by = join_by(date,varid))
  return(out)
})

cpop_summary <- as.data.frame(do.call(rbind, cpop_summary_list))
cpop_summary %>% 
  mutate(outlier = factor(1 + (cpop<0.01), labels = c("Normal","Outlier"))) %>%
  filter(node == 12) %>%
  ggplot(aes(x = date)) + 
  geom_ribbon(aes(ymin = post.q2.5, ymax = post.q97.5), fill = "gray85") +
  geom_point(aes(y = y, col = outlier), shape = 1, size = 0.7) +
  facet_wrap(~varid, ncol = 1) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

cpop_summary %>% 
  group_by(node, varid) %>% 
  summarize(noutliers = sum(cpop<0.01)) %>% 
  print(n = 33)

fixed_summary %>% filter(node==2)
load("set2gaussian.RData")
draws_df <- out$draws_df
q <- 2
pars <- c(paste0("sigma[",1:q,"]"),"A[1,1]","A[2,1]","A[2,2]","Astar[2,1]",paste0("ell[",1:q,"]"),paste0("tau[",1:q,"]"))
library(bayesplot)
color_scheme_set("brewer-Spectral")
mcmc_trace(draws_df,  pars = pars, facet_args = list(ncol = 3)) + facet_text(size = 15)

ggplot(data = mypred, aes(x = date)) + 
  geom_path(aes(y = y), col = "dimgray", linewidth = 0.5) + 
  facet_wrap(~varid, ncol = 1, labeller = label_bquote(y[.(varid)])) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))
