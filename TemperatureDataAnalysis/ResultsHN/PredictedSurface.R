rm(list=ls())
graphics.off()
library(tidyverse)

load("../ApproximateGLGC/TemperatureDataAnalysis/ResultsHN/NNGP_Temps.RData")
yfitted_summary <- yfitted_summary %>% 
  mutate(xcoord = obsCoords[,1]) %>% 
  mutate(ycoord = obsCoords[,2])
pred_summary <- pred_summary %>% 
  mutate(xcoord = prdCoords[,1]) %>% 
  mutate(ycoord = prdCoords[,2])
ypost_summary_nngp <- rbind(yfitted_summary, pred_summary) %>% select(xcoord,ycoord,post.mean,post.sd) %>% gather(key,value,-xcoord,-ycoord) %>% mutate(key = recode(key, `post.mean` = 2, `post.sd` = 5))
rm(yfitted_summary, pred_summary)

###
load("../ApproximateGLGC/TemperatureDataAnalysis/ResultsHN/logNNGP_Temps.RData")
yfitted_summary <- yfitted_summary %>% 
  mutate(xcoord = obsCoords[,1]) %>% 
  mutate(ycoord = obsCoords[,2])
pred_summary <- pred_summary %>% 
  mutate(xcoord = prdCoords[,1]) %>% 
  mutate(ycoord = prdCoords[,2])
ypost_summary_lognngp <- rbind(yfitted_summary, pred_summary) %>% select(xcoord,ycoord,y,post.sd) %>% gather(key,value,-xcoord,-ycoord) %>% mutate(key = recode(key, `y` = 1, `post.sd` = 4))
rm(yfitted_summary, pred_summary)


####
load("../ApproximateGLGC/TemperatureDataAnalysis/ResultsHN/NNHS2_GLGC_Temps.RData")
yfitted_summary <- yfitted_summary %>% 
  mutate(xcoord = obsCoords[,1]) %>% 
  mutate(ycoord = obsCoords[,2])
pred_summary <- pred_summary %>% 
  mutate(xcoord = prdCoords[,1]) %>% 
  mutate(ycoord = prdCoords[,2])
ypost_summary_glgc <- rbind(yfitted_summary, pred_summary) %>% select(xcoord,ycoord,post.mean,post.sd) %>% gather(key,value,-xcoord,-ycoord) %>% mutate(key = recode(key, `post.mean` = 3, `post.sd` = 6))

ypost_summary <- rbind(ypost_summary_nngp, ypost_summary_lognngp, ypost_summary_glgc)

ggp.mean <- ypost_summary %>% filter(key %in% c(1,2,3)) %>% mutate(key = factor(key, labels = c("Observed","GP: Mean", "GLGC: Mean"))) %>% 
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1, barheight = 10)) +
  facet_wrap(~key) +
  theme_void() +
  theme(strip.text = element_text(size=11, vjust = 0),
        legend.title = element_blank(),
        strip.background = element_blank())

ggp.sd <- ypost_summary %>% 
  filter(key %in% c(4,5,6)) %>% 
  mutate(key = factor(key, labels = c("log GP: Standard Deviation", "GP: Standard Deviation", "GLGC: Standard Deviation"))) %>% 
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1, barheight = 10)) +
  facet_wrap(~key) +
  theme_void() +
  theme(strip.text = element_text(size=11, vjust = 0),
        legend.title = element_blank(),
        strip.background = element_blank())
ggp <- gridExtra::grid.arrange(ggp.mean,ggp.sd, nrow = 2)
ggsave(plot = ggp, filename = "./TemperatureDataAnalysis/ResultsHN/summary_temp_data_analysis_mean_sd_surfaces.png", height = 6, width = 10)

