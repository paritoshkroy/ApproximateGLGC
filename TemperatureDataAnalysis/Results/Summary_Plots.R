rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)
library(tune)

load("./TemperatureDataAnalysis/SelectedData/SelectedStatelliteTemps.rda")
str(selected.sat.temps)
obsCoords <- selected.sat.temps %>% filter(!is.na(MaskTemp)) %>% select(Lon,Lat)
prdCoords <- selected.sat.temps %>% filter(is.na(MaskTemp)) %>% select(Lon,Lat)
load("./TemperatureDataAnalysis/ResultsHN/NNHS1_GLGC_Temps.RData")
fixed_summary
pred_summary
min(fit_summary$ess_bulk)
load("./TemperatureDataAnalysis/ResultsHN/NNGP_Temps.RData")
pred_summary
load("./TemperatureDataAnalysis/ResultsHN/logNNGP_Temps.RData")
pred_summary

ypred.summary_gp <- ypred.summary %>% mutate(Model = 1)
yhat_summary_gp <- yhat_summary %>% mutate(Model = 1)
perfor_metric_gp <- perfor_metric
summary_hyper_gp <- data.frame(summary_hyper[,c("mean","sd","2.5%","97.5%")])
summary_hyper_gp <- summary_hyper_gp %>% rownames_to_column()
summary_hyper_gp <- summary_hyper_gp %>% 
  mutate(`Mean(SD)` = paste0(round(mean,3),"(",round(sd,3),")")) %>%
  mutate(`95% CI` = paste0(round(X2.5.,3),",",round(X97.5.,3))) %>%
  mutate(Model = "GP") %>%
  select(Model, rowname,`Mean(SD)`,"95% CI")
xtable::xtable(summary_hyper_gp)

load("~/McGill/ApproximateGLGCRealBigData/predsummaryGLGCMatern32HSHStempData.rda")
ls()
perfor_metric_glgc <- perfor_metric
ypred.summary_glgc <- ypred.summary %>% mutate(Model = 2)
yhat_summary_glgc <- yhat_summary %>% mutate(Model = 2)
summary_hyper_glgc <- data.frame(summary_hyper[,c("mean","sd","2.5%","97.5%")])
summary_hyper_glgc<- summary_hyper_glgc %>% rownames_to_column()
summary_hyper_glgc <- summary_hyper_glgc %>% 
  mutate(`Mean(SD)` = paste0(round(mean,3),"(",round(sd,3),")")) %>%
  mutate(`95% CI` = paste0(round(X2.5.,3),",",round(X97.5.,3))) %>%
  mutate(Model = "GLGC") %>%
  select(Model, rowname,`Mean(SD)`,"95% CI")
xtable::xtable(summary_hyper_glgc)

perform_metric <- rbind(cbind(Model = "GP", perfor_metric_gp), cbind(Model = "GLGC", perfor_metric_glgc))
perform_metric
xtable::xtable(perform_metric, digits = 3)

## Prediction performance
head(ypred.summary_gp)
ypred.obs <- rbind(ypred.summary_gp, yhat_summary_gp) %>% select(xcoord,ycoord,y) %>% 
  rename(post.mean=y) %>% 
  mutate(post.sd=NA) %>%
  mutate(Model = 0)
ypred.gp <- rbind(ypred.summary_gp, yhat_summary_gp) %>% select(xcoord,ycoord,post.mean,post.sd) %>% mutate(Model = 1)
ypred.glgc <- rbind(ypred.summary_glgc, yhat_summary_glgc) %>% select(xcoord,ycoord,post.mean,post.sd)  %>% mutate(Model = 2)
head(ypred.obs)

ypred.df <- rbind(ypred.obs,ypred.gp,ypred.glgc)
ypred.df <- ypred.df %>% mutate(ModelFactor = factor(Model, labels = c("Observed","GP: Mean","GLGC: Mean")))
gpp.mean <- ggplot(ypred.df, aes(x=xcoord,y=ycoord)) + 
  geom_tile(aes(fill = post.mean)) + 
  scale_fill_distiller(
    palette = "Spectral", 
    guide = guide_colourbar(barwidth = 1, barheight = 10)) +
  facet_wrap(~ModelFactor)+
  theme_void() +
  theme(strip.text = element_text(size = 11),
        legend.title = element_blank(),
        strip.background = element_blank())
ypred.df <- ypred.df %>% mutate(ModelFactor = factor(Model, labels = c("","GP: Standard Deviation","GLGC: Standard Deviation")))
gpp.sd <- ggplot(ypred.df, aes(x=xcoord,y=ycoord)) + 
  geom_tile(aes(fill = post.sd)) + 
  facet_wrap(~ModelFactor) +
  scale_fill_distiller(
    palette = "Spectral", 
    guide = guide_colourbar(barwidth = 1, barheight = 10), na.value = "white") +
  theme_void() +
  theme(strip.text = element_text(size=11, vjust = 0),
        legend.title = element_blank(),
        strip.background = element_blank())
gpp <- gridExtra::grid.arrange(gpp.mean,gpp.sd, nrow = 2)
ggsave(plot = gpp, filename = "summary_temp_data_analysis_mean_sd_surfaces.png", height = 6, width = 10)