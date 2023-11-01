rm(list=ls())
graphics.off()
library(fields)
library(Matrix)
library(tidyverse)
library(magrittr)
library(tidybayes)
library(coda)
library(nleqslv)

fpath <- "/home/ParitoshKRoy/git/ApproximateGLGC/"
source(paste0(fpath,"Rutilities/utility_functions.R"))

##########################################################################
## Load and cropped the original data set to get the data for our analysis
##########################################################################
load(paste0(fpath,"./TemperatureDataAnalysis/SelectedData/AllSatelliteTemps.RData"))
head(all.sat.temps)
table(is.na(all.sat.temps$MaskTemp))  # FALSE are training data

xRange <- c(-95.91153,-91.28381) # original
yRange <- c(34.29519,37.06811) # original

cellSize <-  diff(sort(unique(all.sat.temps$Lat)))[1]
cellSize
diff(sort(unique(all.sat.temps$Lon)))[1]
length(sort(unique(all.sat.temps$Lon))) # 500
length(sort(unique(all.sat.temps$Lat))) # 300

xRange <- c(xRange[1]+290*cellSize, xRange[2]-110*cellSize) # cropped (-93.22208 -92.30395)
yRange <- c(yRange[1]+0*cellSize, yRange[2]-199*cellSize) # cropped (34.29519 35.22259)

diff(xRange)/cellSize
diff(yRange)/cellSize

all.sat.temps <- all.sat.temps %>% mutate(idCropped = (Lon>xRange[1] & Lon<xRange[2]) & (Lat>yRange[1] & Lat<yRange[2]))
with(all.sat.temps, table(idCropped))

selected.sat.temps <- all.sat.temps %>% filter(idCropped == 1)
nsite <- nrow(selected.sat.temps)
nsite

##################################################################################
# The object selected.sat.temps.sat.temps is the data to be used in our analysis
# Now select randomly 3000 locations for modeling
##################################################################################
xMidRange <- mean(xRange)
xMidRange
yMidRange <- mean(yRange)
yMidRange

idMissing <- which(is.na(selected.sat.temps$TrueTemp))
length(idMissing)
idNotMissing <- which(!is.na(selected.sat.temps$TrueTemp))
length(idNotMissing)
nsize <- 1000 #3000
set.seed(1)
idSampled <- sample(idNotMissing, size = nsize, replace = FALSE) 
length(idSampled)
set.seed(NULL)
coords <- selected.sat.temps %>% dplyr::select(Lon,Lat) %>% as.matrix() %>% unname()
scaled.coords <- scale(coords, center = c(xMidRange,yMidRange), scale = c(1,1))
xcoord.mean <- attr(scaled.coords, "scaled:center")[1]
xcoord.mean # -92.76301
ycoord.mean <- attr(scaled.coords, "scaled:center")[2]
ycoord.mean # 34.75889
xcoord.sd <- attr(scaled.coords, "scaled:scale")[1]
xcoord.sd
ycoord.sd <- attr(scaled.coords, "scaled:scale")[2]
ycoord.sd
head(scaled.coords)
apply(coords, 2, range)  # range of the spatial domain
apply(scaled.coords, 2, range)

selected.sat.temps$MaskTemp <- selected.sat.temps$TrueTemp
selected.sat.temps$MaskTemp[-idSampled] <- NA
selected.sat.temps <- selected.sat.temps %>% select(Lon,Lat,MaskTemp,TrueTemp)
table(is.na(selected.sat.temps$MaskTemp)) ## FALSE are the data to be used in  model fitting
nsite
save(selected.sat.temps, coords, scaled.coords, file = paste0(fpath,"./TemperatureDataAnalysis/SelectedData/SelectedStatelliteTemps.rda"))
