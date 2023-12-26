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

fname <- list.files(path = ".", pattern = "*\\.RData", full.names = FALSE)

load("NNGP_Temps.RData")
nngp_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord,post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 2, `post.sd` = 6))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)

load("logNNGP_Temps.RData")
lognngp_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord,post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 3, `post.sd` = 7))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)

load("NNHS2_GLGC_Temps.RData")
nnhs2_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord, post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 4, `post.sd` = 8))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)

load("NNHS2_GLGC_Temps.RData")
obs_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  mutate(ysd = NA) %>%
  select(xcoord,ycoord,y,ysd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `y` = 1, `ysd` = 5))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)


ypost_summary <- rbind(nnhs2_pred_dt, lognngp_pred_dt, nngp_pred_dt, obs_dt)

ggp.mean <- ypost_summary %>% 
  filter(key %in% c(1,2,3,4)) %>% 
  mutate(key = factor(key, labels = c("Observed","GP: Mean", "log GP: Mean","GLGC: Mean"))) %>% 
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1, barheight = 10)) +
  facet_wrap(~key, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size=10, vjust = 0),
        legend.title = element_blank(),
        strip.background = element_blank())


ggp.sd <- ypost_summary %>% 
  filter(key %in% c(5,6,7,8)) %>% 
  mutate(key = factor(key, labels = c("", "GP: SD", "log GP: SD", "GLGC: SD"))) %>% 
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1, barheight = 10), 
                       na.value = "white") +
  facet_wrap(~key, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size=9, vjust = 0),
        legend.title = element_blank(),
        strip.background = element_blank())

ggp <- gridExtra::grid.arrange(ggp.mean,ggp.sd, nrow = 2)
ggsave(plot = ggp, filename = "./mean_sd_surfaces_temdata.png", height = 7, width = 11)


