rm(list=ls())
graphics.off()
library(tidyverse)
library(ggplot2)
library(magrittr)
library(coda)
library(fields)
library(lubridate)
library(scoringRules)

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNGP_Temps.RData")
nngp_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord,post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 2, `post.sd` = 6))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/logNNGP_Temps.RData")
lognngp_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord,post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 3, `post.sd` = 7))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)


load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")

nnhs2_pred_dt <- rbind(
  yfitted_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2])) %>%
  select(xcoord,ycoord, post.mean,post.sd) %>% 
  gather(key,value,-xcoord,-ycoord) %>% 
  mutate(key = recode(key, `post.mean` = 4, `post.sd` = 8))
rm(yfitted_summary,pred_summary,obsCoords,prdCoords)

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")
obs_dt <- yfitted_summary %>% 
  mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>%
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
ggsave(plot = ggp, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/mean_sd_surfaces_temdata.png", height = 7, width = 11)

### Scatter plot for SDs
all.equal(nnhs2_pred_dt %>% select(xcoord,ycoord), lognngp_pred_dt %>% select(xcoord,ycoord))

p1 <- inner_join(
  nngp_pred_dt %>% filter(key == 6) %>% select(-key) %>% rename(nngp = value),
  lognngp_pred_dt %>% filter(key == 7) %>% select(-key) %>% rename(lognngp = value), by = c("xcoord","ycoord")) %>%
  ggplot(aes(x = nngp, y = lognngp)) +
  geom_point(size = 0.25) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25) +
  tune::coord_obs_pred() +
  xlab("NNGP") +
  ylab("log NNGP") +
  ggtitle("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
p1


p2 <- inner_join(
  nngp_pred_dt %>% filter(key == 6) %>% select(-key) %>% rename(nngp = value),
  nnhs2_pred_dt %>% filter(key == 8) %>% select(-key) %>% rename(nnhs = value), by = c("xcoord","ycoord")) %>%
  ggplot(aes(x = nngp, y = nnhs)) +
  geom_point(size = 0.25) + 
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25) +
  tune::coord_obs_pred() +
  xlab("NNGP") +
  ylab("NNHS") +
  ggtitle("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
p2


p3 <- inner_join(
  nnhs2_pred_dt %>% filter(key == 8) %>% select(-key) %>% rename(nnhs = value),
  lognngp_pred_dt %>% filter(key == 7) %>% select(-key) %>% rename(lognngp = value), by = c("xcoord","ycoord")) %>%
  ggplot(aes(x = lognngp, y = nnhs)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25) +
  geom_point(size = 0.25) + 
  tune::coord_obs_pred() +
  xlab("log NNGP") +
  ylab("NNHS") +
  ggtitle("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))

ggp <- gridExtra::grid.arrange(p1,p2,p3, nrow = 1)
ggsave(plot = ggp, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/SDPosteriorPredictiveDistribution.png", height = 4, width = 12)

## Distribution of SD at observed data location and predicted data location
rm(list = ls())
load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNGP_Temps.RData")
nngp <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
nngp <- nngp %>% mutate(Key = factor(Key, labels = c("GP: Observed location","GP: Predicted location")))
p1 <- nngp %>% 
  ggplot(aes(x = post.sd)) + 
  geom_histogram(aes(y = after_stat(density)), fill = NA, col = "dimgray") + 
  facet_wrap(~Key, nrow = 2) +
  ylab("Density") +
  xlab("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))


load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/logNNGP_Temps.RData")
lognngp <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
lognngp <- lognngp %>% mutate(Key = factor(Key, labels = c("log GP: Observed location","log GP: Predicted location")))
p2 <- lognngp %>% 
  ggplot(aes(x = post.sd)) + 
  geom_histogram(aes(y = after_stat(density)), fill = NA, col = "dimgray") + 
  facet_wrap(~Key, nrow = 2) +
  ylab("Density") +
  xlab("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")
nnhs <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
nnhs <- nnhs %>% mutate(Key = factor(Key, labels = c("GLGC: Observed location","GLGC: Predicted location")))
p3 <- nnhs %>% 
  ggplot(aes(x = post.sd)) + 
  geom_histogram(aes(y = after_stat(density)), fill = NA, col = "dimgray") + 
  facet_wrap(~Key, nrow = 2) +
  ylab("Density") +
  xlab("SD of the Posterior Predictive Distribution") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
  
p <- gridExtra::grid.arrange(p1,p2,p3, nrow = 1)
ggsave(plot = p, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/DistributionOfSDbyLocationType.png", height = 6, width = 12)

## Is the high value of SD associated with low value of y
rm(list = ls())
load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNGP_Temps.RData")
nngp <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
nngp <- nngp %>% mutate(Key = factor(Key, labels = c("GP: Observed location","GP: Predicted location")))
p1 <- nngp %>% drop_na() %>%
  ggplot(aes(x = y, y = post.sd)) + 
  geom_point(size = 0.25) + 
  ylim(c(0.25,3.5)) +
  facet_wrap(~Key, ncol = 1) +
  xlab("Temerature") +
  ylab("Posterior SD") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
p1

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/logNNGP_Temps.RData")
lognngp <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
lognngp <- lognngp %>% mutate(Key = factor(Key, labels = c("log GP: Observed location","log GP: Predicted location")))
p2 <- lognngp %>% drop_na() %>%
  ggplot(aes(x = y, y = post.sd)) + 
  geom_point(size = 0.25) + 
  ylim(c(0.25,3.5)) +
  facet_wrap(~Key, ncol = 1) +
  xlab("Temerature") +
  ylab("Posterior SD") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
p2

load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")
nnhs <- rbind(yfitted_summary %>% mutate(Key = 1), pred_summary %>% mutate(Key = 2))
nnhs <- nnhs %>% mutate(Key = factor(Key, labels = c("GLGC: Observed location","GLGC: Predicted location")))
p3 <- nnhs %>% drop_na() %>%
  ggplot(aes(x = y, y = post.sd)) + 
  geom_point(size = 0.25) + 
  ylim(c(0.25,3.5)) +
  facet_wrap(~Key, ncol =1) +
  xlab("Temerature") +
  ylab("Posterior SD") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 11))
p3

p <- gridExtra::grid.arrange(p1,p2,p3, nrow = 1)
ggsave(plot = p, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/AssociationBetweenSDandTemperature.png", height = 6, width = 12)

# High SD and Predictive Z
rm(list = ls())
library(tidyverse)
load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")
ls()
z1_summary
hist(pred_summary$post.sd)
nnhs_pred_dt <- pred_summary %>% 
  mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2]) %>%
  select(xcoord,ycoord,post.mean,post.sd)
nnhs_pred_dt %>% filter(post.sd>0.75 & post.sd<0.7501) %>% distinct(xcoord,ycoord)
nnhs_pred_dt %>% filter(post.sd>1.00 & post.sd<1.00025) %>% distinct(xcoord,ycoord)
nnhs_pred_dt %>% filter(post.sd>1.50 & post.sd<1.51) %>% distinct(xcoord,ycoord)
nnhs_pred_dt %>% filter(post.sd>2.50 & post.sd<2.55) %>% distinct(xcoord,ycoord)
nnhs_pred_dt %>% filter(post.sd>3.00)

selectedLocs <- nnhs_pred_dt %>% 
  filter((post.sd>0.75 & post.sd<0.7501) | (post.sd>1.00 & post.sd<1.00025) | (post.sd>1.50 & post.sd<1.51) | (post.sd>2.50 & post.sd < 2.55) | post.sd>3.0) 
selectedLocs <- selectedLocs %>% mutate(SDlevel = cut(post.sd, breaks = c(0.70,0.8,1.25,1.75,2.75,4), labels = 1:5))
selectedLocs <- selectedLocs %>% mutate(SDlevel = factor(SDlevel, labels = c("0.75","1.00","1.51","2.53","3.42")))
selectedLocs <- selectedLocs %>% arrange(SDlevel)
selectedLocs

## Nearest neighbor
nei_info_pred <- FNN::get.knnx(obsCoords, unname(as.matrix(selectedLocs %>% select(xcoord,ycoord)))
, k = nNeighbors); 
str(nei_info_pred)
nnID <- nei_info_pred$nn.index; str(nnID)
nnDist <- nei_info_pred$nn.dist; str(nnDist)

nn_dt1 <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>% .[nnID[1,],] %>% mutate(Site = 0.75)
nn_dt1
nn_dt2 <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>% .[nnID[2,],] %>% mutate(Site = 1.00)
nn_dt2
nn_dt3 <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>% .[nnID[3,],] %>% mutate(Site = 1.51)
nn_dt3
nn_dt4 <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>% .[nnID[4,],] %>% mutate(Site = 2.53)
nn_dt4
nn_dt5 <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2]) %>% .[nnID[5,],] %>% mutate(Site = 3.42)
nn_dt5


p1 <- nnhs_pred_dt %>%
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = post.sd)) +
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1, barheight = 10), 
                       na.value = "white") +
  geom_point(data = selectedLocs, 
             aes(x = xcoord, y = ycoord), shape = 4, size = 1.2) +
  geom_point(data = nn_dt1, aes(x = xcoord, y = ycoord), size = 0.7, shape = 21) +
  geom_point(data = nn_dt2, aes(x = xcoord, y = ycoord), size = 0.7, shape = 21) +
  geom_point(data = nn_dt3, aes(x = xcoord, y = ycoord), size = 0.7, shape = 21) +
  geom_point(data = nn_dt4, aes(x = xcoord, y = ycoord), size = 0.7, shape = 21) +
  geom_point(data = nn_dt5, aes(x = xcoord, y = ycoord), size = 0.7, shape = 21) +
  xlab("") +
  ylab("") +
  theme_void() +
  theme(legend.position = "none")
p1
p2 <- nn_dt1 %>% mutate(rid = row_number()) %>%
  ggplot(aes(x = factor(rid))) +
  ylim(-2.5,2.5) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25,
                width = 0.1) +
  geom_point(aes(y = post.mean), shape = 21) +
  facet_wrap(~Site, labeller = label_bquote("SD of the Posteriror Predictive Distribution:"~.(Site))) +
  xlab("NN sites") +
  ylab(bquote("Posterior mean (95% CI) for"~z[1]*"("*bold(s)*")")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))
p2
p3 <- nn_dt2 %>% mutate(rid = row_number()) %>%
  ggplot(aes(x = factor(rid))) +
  ylim(-2.5,2.5) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25,
                width = 0.1) +
  geom_point(aes(y = post.mean), shape = 21) +
  facet_wrap(~Site, labeller = label_bquote("SD of the Posteriror Predictive Distribution:"~.(Site))) +
  xlab("NN sites") +
  ylab(bquote("Posterior mean (95% CI) for"~z[1]*"("*bold(s)*")")) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))

p4 <- nn_dt3 %>% mutate(rid = row_number()) %>%
  ggplot(aes(x = factor(rid))) +
  ylim(-2.5,2.5) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25,
                width = 0.1) +
  geom_point(aes(y = post.mean), shape = 21) +
  facet_wrap(~Site, labeller = label_bquote("SD of the Posteriror Predictive Distribution:"~.(Site))) +
  xlab("NN sites") +
  ylab(bquote("Posterior mean (95% CI) for"~z[1]*"("*bold(s)*")")) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))


p5 <- nn_dt4 %>% mutate(rid = row_number()) %>%
  ggplot(aes(x = factor(rid))) +
  ylim(-2.5,2.5) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25,
                width = 0.1) +
  geom_point(aes(y = post.mean), shape = 21) +
  facet_wrap(~Site, labeller = label_bquote("SD of the Posteriror Predictive Distribution:"~.(Site))) +
  xlab("NN sites") +
  ylab(bquote("Posterior mean (95% CI) for"~z[1]*"("*bold(s)*")")) +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))


p6 <- nn_dt5 %>% mutate(rid = row_number()) %>%
  ggplot(aes(x = factor(rid))) +
  ylim(-2.5,2.5) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5), linewidth = 0.25,
                width = 0.1) +
  geom_point(aes(y = post.mean), shape = 21) +
  facet_wrap(~Site, labeller = label_bquote("SD of the Posteriror Predictive Distribution:"~.(Site))) +
  xlab("NN sites") +
  ylab(bquote("Posterior mean (95% CI) for"~z[1]*"("*bold(s)*")")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))

p <- gridExtra::grid.arrange(p1,p2,p3, p4, p5, p6, nrow = 2)
ggsave(plot = p, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/DistributionOfz1NNlocations.png", height = 8, width = 12)


# High SD and observed value
rm(list = ls())
library(tidyverse)
load("~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/NNHS3_GLGC_Temps.RData")
ls()
z1_summary
yfitted_summary
nnhs_ydt <- rbind(
  yfitted_summary %>% 
    mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2], Type = 1),
  pred_summary %>% mutate(xcoord = prdCoords[,1], ycoord = prdCoords[,2], Type = 2)) %>%
  select(Type,xcoord,ycoord,post.mean,post.sd,post.q2.5, post.q50, post.q97.5,y)
str(nnhs_ydt)


nnhs_z1dt <- z1_summary %>% mutate(xcoord = obsCoords[,1], ycoord = obsCoords[,2])
str(nnhs_z1dt)

p1 <- nnhs_ydt %>% filter(Type == 1) %>%
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = y)) +
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1.1, barheight = 10), 
                       na.value = "white") +
  geom_point(data = nnhs_ydt %>% filter(Type == 2 & y <= 33),
             aes(x = xcoord, y = ycoord), shape = 4, alpha = 0.5, size = 0.7) +
  geom_point(data = nnhs_ydt %>% filter(Type == 2 & post.sd >= 2),
             aes(x = xcoord, y = ycoord), shape = 21, size = 2) +
  theme_void() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

p2 <- nnhs_z1dt %>% 
  ggplot(aes(x = xcoord, y = ycoord)) + 
  geom_tile(aes(fill = post.mean)) +
  scale_fill_distiller(palette = "Spectral",  
                       guide = guide_colourbar(barwidth = 1.1, barheight = 10), 
                       na.value = "white") +
  geom_point(data = nnhs_ydt %>% filter(Type == 2 & y <= 33),
             aes(x = xcoord, y = ycoord), shape = 4, alpha = 0.5, size = 0.7) +
  geom_point(data = nnhs_ydt %>% filter(Type == 2 & post.sd >= 2),
             aes(x = xcoord, y = ycoord), shape = 21, size = 2) +
  theme_void() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

p3 <- nnhs_ydt %>% 
  filter(Type == 2 & y <= 33) %>%
  mutate(highSD = factor(1+post.sd>2, labels = c("SD <2 ","SD >= 2"))) %>%
  mutate(rid = row_number()) %>%
  ggplot(aes(x = rid, group = highSD)) +
  geom_errorbar(aes(ymin = post.q2.5, ymax = post.q97.5, col = highSD), 
                width = 0.1, linewidth = 0.50) +
  xlab(bquote("Locations with temperature" <= 33)) +
  ylab("Posterior 95% CI for posterior prediction") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.15,0.15),
        panel.grid = element_blank())
p <- gridExtra::grid.arrange(p1,p2,p3,nrow = 1)
ggsave(plot = p, filename = "~/git/ApproximateGLGC/TemperatureDataAnalysis/Round2RevisionPlots/ObseredLocationsLowTemperatureHighPosteriorPredictiveSD.png", height = 4, width = 12)


tmp_dt <- inner_join(nnhs_z1dt %>% rename(z.post.q2.5 = post.q2.5, z.post.q97.5 = post.q97.5) %>% select(xcoord,ycoord,z.post.q2.5,z.post.q97.5),
           nnhs_ydt %>% rename(y.post.q2.5 = post.q2.5, y.post.q97.5 = post.q97.5) %>% select(xcoord,ycoord,y.post.q2.5,y.post.q97.5,post.sd,y))



