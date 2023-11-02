rm(list=ls())
library(tidyverse)
library(sf)
library(lubridate)

fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"

shp <- read_sf(paste0(fpath,"SSTempDataAnalysis/MediterraneanSeaShapeFiles/MediterraneanSea.shp"))
ggplot(shp) + geom_sf(aes(fill = name))
shp$name
Eastern_Mediterranean_Sea <- c("Mediterranean Sea - Eastern Basin", "Aegean Sea", "Ionian Sea" ,"Adriatic Sea")
shp <- shp %>% filter(name %in% Eastern_Mediterranean_Sea)
str(shp)
ggplot(shp) + geom_sf(aes(fill = name))

dt <- read_csv(paste0(fpath,"SSTempDataAnalysis/SelectedData/3504026.csv"))
dt %>% group_by(SST_MM) %>% count()
dt <- dt %>% filter(!is.na(SST_MM))
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
ggplot(shp) + geom_sf(fill = NA) + geom_sf(data = dt_sf, size = 0.1)

dt_shp_sf <- st_filter(dt_sf, shp, .predicate =  st_intersects)
dt_shp_sf <- dt_shp_sf %>% mutate(temp = (temp-32)/1.8) # degree F to degree C
str(dt_shp_sf)
dt_shp_sf %>% distinct(lon,lat)
ggplot(shp) + 
  geom_sf(fill = NA) + 
  geom_sf(data = dt_shp_sf, aes(col = temp), size = 1) +
  scale_color_distiller(palette = "Spectral")
dt_shp_sf %>% ggplot(aes(x = temp)) + geom_density()

eastern_msst <- dt_shp_sf %>% st_drop_geometry()
nsite <- nrow(eastern_msst); nsite
psize <- 100
nsize <- nsite - psize; nsize
set.seed(123)
idSampled <- sample.int(n = nsite, size = nsize, replace = FALSE)
set.seed(NULL)
####################################################################################
## Changing the spatial domain
####################################################################################
coords_Orig <- eastern_msst %>% select(lon,lat) %>% as.matrix() %>% unname()
distMat_Orig <- fields::rdist(coords_Orig)
distVec_Orig <- distMat_Orig[lower.tri(distMat_Orig, diag = FALSE)]
hist(distVec_Orig)

apply(eastern_msst[,c("lon","lat")], 2, range)  # range of the spatial domain
eastern_msst <- eastern_msst %>% 
  mutate(relocateLon = lon - mean(range(lon))) %>%
  mutate(relocateLat = lat - mean(range(lat)))
apply(eastern_msst[,c("relocateLon","relocateLat")], 2, range)

coords <- eastern_msst %>% select(lon,lat) %>% as.matrix() %>% unname()
scaled.coords <- eastern_msst %>% select(relocateLon,relocateLat) %>% as.matrix() %>% unname()
distMat <- fields::rdist(coords)
distVec <- distMat[lower.tri(distMat, diag = FALSE)]

save(eastern_msst, coords, scaled.coords, file = paste0(fpath,"SSTempDataAnalysis/SelectedData/SSTempDataPreparation.rda"))

