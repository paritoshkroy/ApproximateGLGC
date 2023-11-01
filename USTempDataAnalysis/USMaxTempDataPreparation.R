rm(list=ls())
graphics.off()
library(sf)

us_df <- read_csv("./USMaxTempDataAnalysis/2015 USA Weather Data FINAL.csv")
str(us_df)
hist(us_df$`LONGITUDE""`)
library(lubridate)
load("~/McGill/US-EPA-Data/USborder-spatial-polygon.RData")
ls()

class(USborder)
us_border_sfc <- st_as_sfc(USborder)
st_crs(us_border_sfc)
#us_border_sfc <- st_transform(us_border_sfc, crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
ggplot(us_border_sfc) + geom_sf(fill = NA, col = 1)

us_temp <- read_csv("~/Downloads/daily_TEMP_2018.zip")
us_temp <- us_temp %>% filter(`Parameter Code` == 62101)
unique(us_temp$`Parameter Name`)
unique(us_temp$`Units of Measure`)
us_temp <- us_temp %>% mutate(date = as_date(`Date Local`)) %>%
  mutate(month = month(date)) %>%
  mutate(day = day(date)) %>% 
  mutate(year = year(date)) %>%
  rename(outdoor_tem_degF = `Arithmetic Mean`) %>%
  select(Latitude,Longitude,date,year,month,day,outdoor_tem_degF)
us_temp
us_summer_temp <- us_temp %>% filter(month %in% c(9,10,11))
us_summer_temp
us_summer_tmax <- us_summer_temp %>%
  group_by(Latitude, Longitude) %>%
  summarise(tmax = max(outdoor_tem_degF)) %>%
  ungroup()
us_summer_tmax

ggplot(us_summer_tmax, aes(x = tmax)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 31, 
                 col = "dimgray", fill = NA) + 
  geom_density()
boxplot(us_summer_tmax$tmax)

##?????????????? Need to check 
us_tmax_raw <- us_temp
us_tmax_sf_raw <- st_as_sf(us_tmax_raw, coords = c("Longitude","Latitude"), crs = st_crs(us_border_sfc))

id <- st_intersects(us_border_sfc, us_tmax_sf_raw)
str(id)
ggplot(us_tmax_sf_raw[id[[1]],]) + geom_sf()
us_tmax_sf_within <- us_tmax_sf_raw[id[[1]],]

ggplot() + geom_sf(data = us_border_sfc, fill = NA, col = 2) +
  geom_sf(data = us_tmax_sf_within, aes(col = tmax)) + 
  scale_color_viridis() +
  theme_void() +
  theme(legend.title = element_blank())

library(elevatr)
us_tmax_sf <- get_elev_point(us_tmax_sf_within, prj = st_crs(us_tmax_sf_within), src = "epqs")
str(us_tmax_sf)

ggplot() + geom_sf(data = us_border_sfc, fill = NA, col = 2) +
  geom_sf(data = us_tmax_sf, aes(col = elevation)) + 
  scale_color_viridis() +
  theme_void() +
  theme(legend.title = element_blank())

save(us_tmax_sf, us_border_sfc, file = "us_tamx_mar2018.rda")


