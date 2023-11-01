rm(list=ls())
graphics.off()
library(tidyverse)
library(sp)
data(meuse)
str(meuse)
library(sf)
meuse.sf <- data.frame(meuse) %>% 
  mutate(x=x/1000) %>%
  mutate(y=y/1000) %>%
  st_as_sf(coords=c("x","y"))
st_crs(meuse.sf) <- 28992
meuse.border <- st_convex_hull(st_union(st_buffer(meuse.sf, dist = 0.15)))
ggplot(meuse.border) + geom_sf(fill = NA) +
  geom_sf(data = meuse.sf, aes(size = zinc)) + 
  theme_bw() + 
  theme(panel.grid = element_blank())

meuse.grid <- meuse.border %>%
  st_make_grid(cellsize = 50/1000, what = "centers") %>% # grid of points
  st_intersection(meuse.border) %>% # only within the polygon
  st_as_sf()
meuse.grid %>% dim()
data(meuse.grid) # ignoring own grid location and used the given grid location in sp package
meuse.grid.sf <- data.frame(meuse.grid) %>% 
  mutate(x=x/1000) %>%
  mutate(y=y/1000) %>%
  st_as_sf(coords=c("x","y"))
st_crs(meuse.grid.sf) <- 28992

data(meuse.riv)
meuse.riv.sf <- data.frame(meuse.riv) %>%  
  mutate(x=X1/1000) %>%
  mutate(y=X2/1000) %>%
  st_as_sf(coords=c("x","y"))
st_crs(meuse.riv.sf) <- 28992

ggplot(meuse.grid.sf) + 
  geom_sf(color = "dimgray", shape=3, size=1) + 
  geom_sf(data = meuse.sf, aes(size=zinc/100), colour = "blue", shape = 1) +
  theme_bw()+
  theme(panel.grid = element_blank())


obsCoords <- st_coordinates(meuse.sf) %>% as.matrix() %>% .[,1:2]
colnames(obsCoords) <- c("longitude","latitude")
dim(obsCoords)
head(obsCoords)

coords0 <- st_coordinates(meuse.grid.sf) %>% as.matrix() %>% .[,1:2]
colnames(coords0) <- c("longitude","latitude")
dim(coords0)
head(coords0)

# EDA
# because in the prediction location we have 
# the distance to river in a scaled version within [0,1]
meuse.sf <- meuse.sf %>% mutate(y=zinc/100)
meuse.sf <- meuse.sf %>% mutate(dist=dist.m/max(dist.m))
meuse.sf <- meuse.sf %>% mutate(tdist=sqrt(dist))
meuse.sf <- meuse.sf %>% mutate(fdummy=if_else(ffreq==1,1,0))

hist(log(meuse.sf$y))

ggplot(meuse.sf) + 
  geom_boxplot(aes(x=factor(fdummy),y=y)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# extracting the data for modeling
y <- meuse.sf$y
tdist <- meuse.sf$tdist
fdummy <- meuse.sf$fdummy
X <- cbind(tdist,fdummy)
head(X)
p <- ncol(X)
nsize <- length(y); nsize
