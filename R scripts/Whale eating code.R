
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

bsw_whales_1_breed_season = bsw_whales_1_consumption_day * 30
bsw_whales_2_breed_season = bsw_whales_2_consumption_day * 30
bsw_whales_3_breed_season = bsw_whales_3_consumption_day * 30

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

bse_whales_1_breed_season = bse_whales_1_consumption_day * 30
bse_whales_2_breed_season = bse_whales_2_consumption_day * 30
bse_whales_3_breed_season = bse_whales_3_consumption_day * 30

bse_whales_total_breed_season = bse_whales_1_breed_season + bse_whales_2_breed_season + bse_whales_3_breed_season


TOTAL_BS_WHALE = bse_whales_total_breed_season + bsw_whales_total_breed_season

####Now we read in the catch data from Watters et al, subset the two SSMU, sum the catch for the same periods, and compare----
catch = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APBSE", "APBSW"))%>%#, CalendarYear >= "2010") %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

####back of the fag packt reverse calc for male AFS to hit the 155kt 48.1 trigger----
#uses Boyd et al. consumption estimate of a 4yr old male eating 3tonnes/year, or an average of 8kg/day

quota = 31000
seal_consumption = 0.008*90 #8kg per day for 90days, or between January 1st and March 31st each year
number_seals = (quota/seal_consumption)
number_seals


catch_and_seals = catch %>%
  dplyr::mutate(equivalent_AFS_n = TotalCatch/0.72)