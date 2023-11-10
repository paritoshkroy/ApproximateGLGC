rm(list=ls())
library(tidyverse)
library(sf)
library(lubridate)
fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"

shp <- read_sf(paste0(fpath,"SSTempDataAnalysis/MediterraneanSeaShapeFiles/MediterraneanSea.shp"))
EastMedSea <- c("Mediterranean Sea - Eastern Basin", "Aegean Sea", "Ionian Sea" ,"Adriatic Sea")
shp <- shp %>% filter(name %in% EastMedSea)
ggplot(shp) + geom_sf(aes(fill = name))
shp$name
ggplot(shp) + geom_sf(aes(fill = name)) + theme_void() + theme(legend.title = element_blank())
st_bbox(shp)

dt <- read_csv(paste0(fpath,"SSTempDataAnalysis/SelectedData/3504026.csv"))
dt %>% group_by(SST_MM) %>% count()
dt <- dt %>% filter(!is.na(SST_MM)) %>% filter(SST_MM != 0)
dt <- dt %>% mutate(date = as.Date(DATE)) %>% select(date,LONGITUDE,LATITUDE,SEA_SURF_TEMP)
dt <- dt %>% drop_na()
dt %>% distinct(LONGITUDE,LATITUDE)
dt %>% select(LONGITUDE,LATITUDE) %>% duplicated()
dt %>% arrange(date,LONGITUDE,LATITUDE)
dt %>% distinct(date) %>% arrange(date)
dt %>% group_by(date) %>% count()
dt %>% group_by(LONGITUDE,LATITUDE) %>% count()
dt <- dt %>% group_by(LONGITUDE,LATITUDE) %>%
  summarise(temp = median(SEA_SURF_TEMP)) %>%
  ungroup()
str(dt)
dt %>% distinct(LONGITUDE,LATITUDE)
dt <- dt %>% mutate(lon = LONGITUDE,lat = LATITUDE)
dt <- dt %>% filter(!(lon>32 & lat <30.5))
dt_sf <- st_as_sf(dt, coords = c("LONGITUDE","LATITUDE"))
st_crs(dt_sf) <- st_crs(shp)
st_bbox(shp)
ggplot(shp) + geom_sf(fill = NA) + geom_sf(data = dt_sf, size = 0.25)

dt_shp_sf <- st_filter(dt_sf, shp, .predicate =  st_intersects)
dt_shp_sf <- dt_shp_sf %>% mutate(temp = (temp-32)/1.8) # degree F to degree C
str(dt_shp_sf)
dt_shp_sf %>% distinct(lon,lat)
ggplot(shp) + 
  geom_sf(fill = NA) + 
  geom_sf(data = dt_shp_sf, aes(col = temp), shape = 20, size = 1) +
  scale_color_distiller(palette = "Spectral") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key.height = unit(1, 'cm'), 
        legend.position = "right")
ggsave(filename = "./SSTempDataAnalysis/SelectedData/SSTempDataLocations.png", height = 4, width = 6)

dt_shp_sf <- dt_shp_sf %>% 
  mutate(resid = as.numeric(residuals(lm(temp~lon+lat)))) 
dt_shp_sf %>%
  ggplot(aes(x = resid)) +
  geom_histogram(aes(y = after_stat(density)), fill = NA, col = "dimgray", bins = 21) +
  geom_density() +
  xlab("Residuals") +
  ylab("Density") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid = element_blank())
ggsave(filename = "./SSTempDataAnalysis/SelectedData/SSTempResidualsDistribution.png", height = 4, width = 6)

msst_df <- dt_shp_sf %>% st_drop_geometry()
nsite <- nrow(msst_df); nsite
nsize <- 500 #ceiling(nsite*0.50); nsize #1500
psize <- nsite - nsize; psize

set.seed(123)
idSampled <- sample.int(n = nsite, size = nsize, replace = FALSE)
set.seed(NULL)
####################################################################################
## Changing the spatial domain
####################################################################################
coords_Orig <- msst_df %>% select(lon,lat) %>% as.matrix() %>% unname()
distMat_Orig <- fields::rdist(coords_Orig)
distVec_Orig <- distMat_Orig[lower.tri(distMat_Orig, diag = FALSE)]
hist(distVec_Orig)

apply(msst_df[,c("lon","lat")], 2, range)  # range of the spatial domain
msst_df <- msst_df %>% 
  mutate(relocateLon = lon - mean(range(lon))) %>%
  mutate(relocateLat = lat - mean(range(lat)))

apply(msst_df[,c("relocateLon","relocateLat")], 2, range)

coords <- msst_df %>% select(lon,lat) %>% as.matrix() %>% unname()
scaled.coords <- msst_df %>% select(relocateLon,relocateLat) %>% as.matrix() %>% unname()
distMat <- fields::rdist(coords)
distVec <- distMat[lower.tri(distMat, diag = FALSE)]
abline(v = quantile(distVec, probs = seq(0,1,l = 21)), col = 2)
quantile(distVec, probs = seq(0,1,l=21))

## For minimum m1 and m2 for the HSGP
Lstar <- as.vector(apply(apply(msst_df[,c("relocateLon","relocateLat")], 2, range),2,max)); Lstar
quantile(distVec, probs = c(1,2.5,5)/100)
minimum_identifiable_lscale <- 1.20; minimum_identifiable_lscale
c <- max(round(1.2 + minimum_identifiable_lscale/Lstar,digits = 1)); c
m1 <- ceiling(3.42 * c/(minimum_identifiable_lscale/Lstar[1])); m1
m2 <- ceiling(3.42 * c/(minimum_identifiable_lscale/Lstar[2])); m2
m1*m2

save(m1, m2, idSampled, msst_df, coords, scaled.coords, file = paste0(fpath,"SSTempDataAnalysis/SelectedData/SSTempDataPreparation.rda"))
