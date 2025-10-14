
library(sf)
library(terra)
library(tidyverse)


whales_1 = terra::rast("/Volumes/GodDrive/Quantarctica3/IAATO DSHM/Cruise 1.grd")
whales_2 = terra::rast("/Volumes/GodDrive/Quantarctica3/IAATO DSHM/Cruise 2.grd")
whales_3 = terra::rast("/Volumes/GodDrive/Quantarctica3/IAATO DSHM/Cruise 3.grd")

whales_1 <- ifel(is.na(whales_1), 0, whales_1)
whales_2 <- ifel(is.na(whales_2), 0, whales_2)
whales_3 <- ifel(is.na(whales_3), 0, whales_3)

ssmu = sf::read_sf("/Volumes/GodDrive/Quantarctica3/EnvironmentalManagement/CCAMLR/CCAMLR_SSMUs.shp") %>%
  sf::st_transform(crs = sf::st_crs(whales_1))
  
bsw = ssmu[5, ] %>%
  terra::vect() 
bse = ssmu[6,] %>%
  terra::vect()
####BS West----
bsw_whales_1 = terra::extract(whales_1, bsw, fun = sum)
bsw_whales_2 = terra::extract(whales_2, bsw, fun = sum)
bsw_whales_3 = terra::extract(whales_3, bsw, fun = sum)

bsw_whales_1
bsw_whales_2
bsw_whales_3

bsw_whales_1_consumption_day = bsw_whales_1$z * 0.632
bsw_whales_2_consumption_day = bsw_whales_2$z * 0.632
bsw_whales_3_consumption_day = bsw_whales_3$z * 0.632

bsw_whales_1_breed_season = bsw_whales_1_consumption_day * 90
bsw_whales_2_breed_season = bsw_whales_2_consumption_day * 90
bsw_whales_3_breed_season = bsw_whales_3_consumption_day * 90

bsw_whales_total_breed_season = bsw_whales_1_breed_season + bsw_whales_2_breed_season + bsw_whales_3_breed_season
####BS EAST----
bse_whales_1 = terra::extract(whales_1, bse, fun = sum)
bse_whales_2 = terra::extract(whales_2, bse, fun = sum)
bse_whales_3 = terra::extract(whales_3, bse, fun = sum)

bse_whales_1
bse_whales_2
bse_whales_3

bse_whales_1_consumption_day = bse_whales_1$z * 0.632
bse_whales_2_consumption_day = bse_whales_2$z * 0.632
bse_whales_3_consumption_day = bse_whales_3$z * 0.632

bse_whales_1_breed_season = bse_whales_1_consumption_day * 90
bse_whales_2_breed_season = bse_whales_2_consumption_day * 90
bse_whales_3_breed_season = bse_whales_3_consumption_day * 90

bse_whales_total_breed_season = bse_whales_1_breed_season + bse_whales_2_breed_season + bse_whales_3_breed_season


TOTAL_BS_WHALE = bse_whales_total_breed_season + bsw_whales_total_breed_season
