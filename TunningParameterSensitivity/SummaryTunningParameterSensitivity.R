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
load("./TunningParameterSensitivity/TunningParameterSensitivitySummary.RData")
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

fixed_summary <- fixed_summary %>% separate(Method, into = c("x","y"), sep = "_") %>% 
  mutate(x = recode(x, Full="Full0")) %>% 
  mutate(y = recode(y, GLGC="GLGC0")) %>%
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,LS,everything()) %>%
  mutate(ellfactor = factor(LS, labels = seq(0.1,1,l=10))) %>%
  mutate(m1factor = factor(m1, labels = c("Full",sort(unique(m1))[-1]))) %>%
  mutate(m2factor = factor(m2, labels = c("Full",sort(unique(m2))[-1])))

fixed_summary %>% 
  filter(variable %in% "ell2") %>% 
  ggplot(aes(x = ellfactor, group = m1factor)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = m1factor), 
                position=position_dodge(0.5), stat="identity", 
                width = 0.25) +
  facet_wrap(~Method, nrow = 2) +
  labs(x = bquote("\u2113"[2])) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))


fixed_summary %>% 
  filter(variable %in% "ell2") %>% 
  ggplot(aes(x = factor(LS), group = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`, color = Method), 
                position=position_dodge(0.5), stat="identity", 
                width = 0.25) +
  #facet_wrap(~LS, nrow = 2, labeller = label_bquote("\u2113"[2]~"="~.(LS))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 12))

fixed_summary %>% 
  filter(variable %in% "gamma") %>% 
  filter(LS == 7) %>%
  ggplot(aes(x = Method)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25) +
  geom_point(aes(y = `50%`)) +
  geom_hline(yintercept = 1.5) +
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
scores_df <- scores_df %>% mutate(x = recode(x, Full="Full0")) %>% mutate(y = recode(y, GLGC="GLGC0"))
scores_df <- scores_df %>% 
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,LS,everything())
scores_df <- scores_df %>% mutate(ell = factor(LS, labels = seq(0.1,1,l=10)))
scores_df <- scores_df %>% mutate(m1factor = factor(m1, labels = c("Full",sort(unique(m1))[-1])))
scores_df <- scores_df %>% mutate(m2factor = factor(m2, labels = c("Full",sort(unique(m2))[-1])))
scores_df %>% select(m1factor)

#scores_df <- scores_df  %>% filter(LS <= 4)
## Energy Score
library(gridExtra)
maxmin_escore <- range(scores_df$ES)
escore_plots <- grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN")) %>% 
    ggplot(aes(x = ell, y = ES)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNNN") +
    xlab("\u2113")+
    ylab("Energy Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_escore)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS")) %>% 
    ggplot(aes(x = ell, y = ES)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNHS") +
    xlab("\u2113")+
    ylab("Energy Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_escore)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","HSHS")) %>% 
    ggplot(aes(x = ell, y = ES)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("HSHS") +
    xlab("\u2113")+
    ylab("Energy Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_escore)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)), 
  
  ncol = 3)

## log Score
maxmin_logs <- range(scores_df$logs)
logs_plots <- grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN")) %>% 
    ggplot(aes(x = ell, y = logs)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNNN") +
    xlab("\u2113")+
    ylab("log Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_logs)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS")) %>% 
    ggplot(aes(x = ell, y = logs)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNHS") +
    xlab("\u2113")+
    ylab("log Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_logs)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","HSHS")) %>% 
    ggplot(aes(x = ell, y = logs)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("HSHS") +
    xlab("\u2113")+
    ylab("log Score")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_logs)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)), 
  
  ncol = 3)

## Interval Score
maxmin_is <- range(scores_df$IS)
iscore_plots <- grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN")) %>% 
    ggplot(aes(x = ell, y =IS)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNNN") +
    xlab("\u2113")+
    ylab("IS")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_is)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS")) %>% 
    ggplot(aes(x = ell, y = IS)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("NNHS") +
    xlab("\u2113")+
    ylab("IS")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_is)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","HSHS")) %>% 
    ggplot(aes(x = ell, y = IS)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 2) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.7) +
    ggtitle("HSHS") +
    xlab("\u2113")+
    ylab("IS")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    ylim(maxmin_is)+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 12), 
          legend.position = c(0.8,0.6), 
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 11),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12)), 
  
  ncol = 3)

multiscores <- grid.arrange(iscore_plots, escore_plots, logs_plots, ncol = 1)
ggsave(multiscores, filename = "multiscores.png", height = 5, width = 10)

names(scores_df)
multiscores <- grid.arrange(
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNNN")) %>% 
    select(ell,m1factor,m2factor,MAE,RMSE,IS,CRPS,CVG) %>%
    gather(key,value,-c(ell,m1factor,m2factor)) %>%
    ggplot(aes(x = ell, y = value)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 0.35) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.35) +
    facet_wrap(~key, scales = "free_y", nrow = 1) +
    xlab("")+
    ylab("NNNN")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 10), 
          legend.position = c(0.95,0.6), 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.25, "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 6),
          axis.ticks.x = element_blank()),
  
  scores_df %>% 
    filter(Method %in% c("FullGLGC","NNHS")) %>% 
    select(ell,m1factor,m2factor,MAE,RMSE,IS,CRPS,CVG) %>%
    gather(key,value,-c(ell,m1factor,m2factor)) %>%
    ggplot(aes(x = ell, y = value)) + 
    geom_point(aes(color= m1factor), shape = 1, size = 0.35) + 
    geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.35) +
    facet_wrap(~key, scales = "free_y", nrow = 1) +
    xlab("")+
    ylab("NNHS")+
    labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          legend.position = c(0.95,0.6), 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.25, "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 6),
          axis.ticks.x = element_blank()),
  
  scores_df %>% 
  filter(Method %in% c("FullGLGC","HSHS")) %>% 
  select(ell,m1factor,m2factor,MAE,RMSE,IS,CRPS,CVG) %>%
  gather(key,value,-c(ell,m1factor,m2factor)) %>%
  ggplot(aes(x = ell, y = value)) + 
  geom_point(aes(color= m1factor), shape = 1, size = 0.35) + 
  geom_line(aes(color = m1factor, group = m1factor), linewidth = 0.35) +
  facet_wrap(~key, scales = "free_y", nrow = 1) +
  xlab("\u2113")+
  ylab("HSHS")+
  labs(color = expression(m*"/"*m[1]*"/"*m[2]))+
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = c(0.95,0.6), 
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm"),
        axis.title.x = element_text(size = 10),
        axis.text = element_text(size = 6)),
  
  nrow = 3)
ggsave(multiscores, filename = "multiscores.png", height = 6, width = 10)

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

save(scores_list, fixed_summary_list, file = "TunningParameterSensitivitySummary.RData")