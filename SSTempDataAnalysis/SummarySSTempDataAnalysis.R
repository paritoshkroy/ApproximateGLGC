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
fname <- list.files(pattern = "*\\.RData", full.names = FALSE)

fixed_summary_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  fixed_summary <- fixed_summary %>% 
    mutate(Model = fname[node]) %>% 
    select(Model,everything())
  return(fixed_summary)
})
fixed_summary <- do.call(rbind, fixed_summary_list)
fixed_summary <- fixed_summary %>% separate(Model, into = c("x","y"), sep = "_G") %>%
  mutate(x = recode(x, Full="Full0_Full0")) %>% select(-y) %>% rename(Model=x) %>%
  separate(Model, into = c("x","y"), sep = "_") %>%
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  mutate(Method = recode(Method, FullFull="Full")) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,everything())

fixed_summary <- fixed_summary %>% mutate(MethodCharc = paste0(Method,"(",m1,",",m2,")")) 
fixed_summary <- fixed_summary %>% 
  mutate(MethodFactor = recode(MethodCharc, "Full(0,0)" = 1, "NNNN(15,15)" = 2, "NNNN(10,10)" = 3, "NNNN(6,6)" = 4, "HSHS(51,28)" = 5, "HSHS(44,24)" = 6)) 
fixed_summary <- fixed_summary %>% 
  mutate(MethodFactor = factor(MethodFactor, labels = c("Full", "NNNN(15,15)", "NNNN(10,10)", "NNNN(6,6)", "HSHS(51,28)", "HSHS(44,24)")))
fixed_summary <- fixed_summary %>% 
  mutate(Pars = recode(variable, `theta[1]`=1, `theta[2]` = 2, `theta[3]` = 3, gamma = 4, sigma1 = 5, sigma2 = 6, ell1 = 7, ell2 = 8, tau = 9)) %>%
  mutate(Pars = factor(Pars, labels = c("theta[1]","theta[2]","theta[3]","gamma","sigma[1]","sigma[2]","\u2113[1]","\u2113[2]","tau"))) 

ggplot(fixed_summary, aes(x = MethodFactor)) + 
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.25,
                position=position_dodge(0.5), stat="identity") + 
  facet_wrap(~Pars, scales = "free_y", labeller = label_parsed, nrow = 3) +
  theme_bw() +
  xlab("") +
  ylab("Posterior median (95% CI)") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = "FixedParsSummaries.png", height = 5, width = 9)

fixed_summary$MethodFactor

###################################################################
## Predictive Scores
###################################################################
scores_list <- lapply(1:length(fname), function(node){
  load(fname[node])
  if(exists("scores_df")) return(scores_df %>% mutate(Model = fname[node]))
})
scores_df <- do.call(rbind, scores_list)
scores_df <- scores_df %>% separate(Model, into = c("x","y"), sep = "_G") %>%
  mutate(x = recode(x, Full="Full0_Full0")) %>% select(-y) %>% rename(Model=x) %>%
  separate(Model, into = c("x","y"), sep = "_") %>%
  mutate(C1 = gsub("[^A-Za-z]","",x)) %>%
  mutate(C2 = gsub("[^A-Za-z]","",y)) %>%
  mutate(m1 = as.numeric(gsub(".*?([0-9]+).*", "\\1", x))) %>%
  mutate(m2 = as.numeric(gsub(".*?([0-9]+).*", "\\1", y))) %>%
  mutate(Method = paste0(C1,C2)) %>%
  mutate(Method = recode(Method, FullFull="Full")) %>%
  select(-x,-y) %>%
  select(Method,m1,m2,everything())

scores_df <- scores_df %>% mutate(MethodCharc = paste0(Method,"(",m1,",",m2,")")) 
scores_df <- scores_df %>% 
  mutate(MethodFactor = recode(MethodCharc, "Full(0,0)" = 1, "NNNN(15,15)" = 2, "NNNN(10,10)" = 3, "NNNN(6,6)" = 4, "HSHS(51,28)" = 5, "HSHS(44,24)" = 6)) 
scores_df <- scores_df %>% 
  mutate(MethodFactor = factor(MethodFactor, labels = c("Full", "NNNN(15,15)", "NNNN(10,10)", "NNNN(6,6)", "HSHS(51,28)", "HSHS(44,24)")))

scores_df %>% select(MethodFactor, MAE, RMSE, CVG, CRPS, IS, ES, logs, `Elapsed Time`) %>% arrange(MethodFactor) %>% xtable::xtable()

##
rm(list=ls())
load("~/git/ApproximateGLGC/SSTempDataAnalysis/NNGP_SST.RData")
str(prdGrid)
NNGPPred <- st_sf(prdGrid) %>% 
  mutate(x = st_coordinates(prdGrid)[,1]) %>% 
  mutate(y = st_coordinates(prdGrid)[,2]) %>% 
  mutate(post_mean = pred_summary$post.mean) %>%
  mutate(post_sd = pred_summary$post.sd) %>%
  mutate(Model = 1)
lpNNGP <- summary(draws_df$lp__)

load("~/git/ApproximateGLGC/SSTempDataAnalysis/logNNGP_SST.RData")
str(prdGrid)
logNNGPPred <- st_sf(prdGrid) %>% 
  mutate(x = st_coordinates(prdGrid)[,1]) %>% 
  mutate(y = st_coordinates(prdGrid)[,2]) %>% 
  mutate(post_mean = pred_summary$post.mean) %>%
  mutate(post_sd = pred_summary$post.sd) %>%
  mutate(Model = 2)
lplogNNGP <- summary(draws_df$lp__)

load("~/git/ApproximateGLGC/SSTempDataAnalysis/HSHS_GLGC_SST.RData")
str(prdGrid)
HSHSPred <- st_sf(prdGrid) %>% 
  mutate(x = st_coordinates(prdGrid)[,1]) %>% 
  mutate(y = st_coordinates(prdGrid)[,2]) %>% 
  mutate(post_mean = pred_summary$post.mean) %>%
  mutate(post_sd = pred_summary$post.sd) %>%
  mutate(Model = 3)
lpHSHS <- summary(draws_df$lp__)

rbind(lpNNGP, lplogNNGP, lpHSHS)

PredDF <- rbind(NNGPPred, logNNGPPred, HSHSPred)
PredDF <- PredDF %>% mutate(Model = factor(Model, labels = c("NNGP","logNNGP","HSHS")))
ggp <- gridExtra::grid.arrange(
  
  ggplot(PredDF) + 
    geom_raster(aes(x = x, y = y, fill = post_mean)) + 
    scale_fill_distiller(palette = 'Spectral') +
    facet_wrap(~Model, nrow = 1) +
    coord_equal()+
    theme_void() +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          legend.key.height = unit(0.7, 'cm'), 
          legend.position = "right"),
  
  ggplot(PredDF) + 
    geom_raster(aes(x = x, y = y, fill = post_sd)) + 
    scale_fill_distiller(palette = 'Spectral') +
    facet_wrap(~Model, nrow = 1) +
    coord_equal()+
    theme_void() +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          legend.key.height = unit(0.7, 'cm'), 
          legend.position = "right"),
  ncol = 1)
ggp
ggsave(plot = ggp, filename = "./SSTempDataAnalysis/SST_Prediction_MeasnSD.png", height = 4, width = 10)


## HSGP
A tibble: 6 × 9
variable    mean     sd   `50%`   `2.5%` `97.5%`  rhat ess_bulk ess_tail
<chr>      <num>  <num>   <num>    <num>   <num> <num>    <num>    <num>
1 sigma     0.506  0.175   0.490   0.203     0.902  1.00     539.     895.
2 tau       1.81   0.0488  1.81    1.71      1.91   1.00    3298.    1543.
3 ell       2.74   0.890   2.60    1.48      4.86   1.00    1945.    1603.
4 theta[1] 18.3    0.249  18.3    17.9      18.8    1.00    2388.    1473.
5 theta[2]  0.0655 0.0382  0.0617  0.00154   0.151  1.00    1638.    1222.
6 theta[3] -0.397  0.0946 -0.405  -0.565    -0.188  1.00    1906.    1401.
