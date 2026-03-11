library(sf)
library(terra)
library(tidyverse)
library(trip)
library(adehabitatHR)

d = read_csv("./Hinke_telemetry/Data/satellite telemetry.csv")

d$Date = as.Date(d$Date, format = "%m/%d/%Y")
d$Year = format(d$Date, "%Y")
d$Month = format(d$Date, "%m")
d$month_day <- format(d$Date, "%m-%d")
d$datetime <- as.POSIXct(paste(d$Date, d$Time), format="%Y-%m-%d %H:%M:%S")

d = d %>%
  dplyr::filter(Loc.Qual %in% c("3", "2", "1", "0", "A"))#, "B"))

#CHINSTRAPS-
####Cape Shirreff----
d_CHPE_CS = d %>%
  dplyr::filter(Spp == "CHPE" & Site == "CS" & month_day < "03-10" & Deployment != "216") %>% #REMOVE the single deployment responsible for the long trip into the Drake
  #dplyr::filter(month_day >= "11-15" | month_day <= "03-10") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(3031)

plot(st_geometry(d_CHPE_CS))

d_CHPE_CS_sp = as(d_CHPE_CS, "Spatial")

d_CHPE_CS_UD = adehabitatHR::kernelUD(d_CHPE_CS_sp, 
                                  h = "href", 
                                  grid = 500, 
                                  kern = "bivnorm", 
                                  same4all = TRUE)
d_CHPE_CS_UD_vert = getverticeshr(d_CHPE_CS_UD, percent = 95)
plot(d_CHPE_CS_UD_vert)

d_CHPE_CS_UD_area = kernel.area(d_CHPE_CS_UD, 
                                percent = c(90),
                                unin = c("m"),
                                unout = c("km2"), standardize = FALSE)

plot(getverticeshr(d_CHPE_CS_UD, percent = 95))

sqnm = d_CHPE_CS_UD_area * 0.29155334959812 
sqnm

####Copacabana----
d_CHPE_COPA = d %>%
  dplyr::filter(Spp == "CHPE" & Site == "COPA" & month_day < "03-10" & Deployment != "216") %>% #REMOVE the single deployment responsible for the long trip into the Drake
#dplyr::filter(month_day >= "11-15" | month_day <= "03-10") %>%
sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(3031)

plot(st_geometry(d_CHPE_COPA))

d_CHPE_COPA_sp = as(d_CHPE_COPA, "Spatial")

d_CHPE_COPA_UD = adehabitatHR::kernelUD(d_CHPE_COPA_sp, 
                                      h = "href", 
                                      grid = 500, 
                                      kern = "bivnorm", 
                                      same4all = TRUE)
d_CHPE_COPA_UD_vert = getverticeshr(d_CHPE_COPA_UD, percent = 95)
plot(d_CHPE_COPA_UD_vert)

d_CHPE_COPA_UD_area = kernel.area(d_CHPE_COPA_UD, 
                                percent = c(90),
                                unin = c("m"),
                                unout = c("km2"), standardize = FALSE)

plot(getverticeshr(d_CHPE_COPA_UD, percent = 95))

sqnm = d_CHPE_COPA_UD_area * 0.29155334959812 
sqnm

#GENTOOS NEXT
####Copacabana----
d_GEPE_COPA = d %>%
  dplyr::filter(Spp == "GEPE" & Site == "COPA") %>% #REMOVE the single deployment responsible for the long trip into the Drake
  #dplyr::filter(month_day >= "11-15" | month_day <= "03-10") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(3031)

plot(st_geometry(d_GEPE_COPA))

d_GEPE_COPA_sp = as(d_GEPE_COPA, "Spatial")

d_GEPE_COPA_UD = adehabitatHR::kernelUD(d_GEPE_COPA_sp, 
                                        h = "href", 
                                        grid = 500, 
                                        kern = "bivnorm", 
                                        same4all = TRUE)
d_GEPE_COPA_UD_vert = getverticeshr(d_GEPE_COPA_UD, percent = 95)
plot(d_GEPE_COPA_UD_vert)

d_GEPE_COPA_UD_area = kernel.area(d_GEPE_COPA_UD, 
                                  percent = c(90),
                                  unin = c("m"),
                                  unout = c("km2"), standardize = FALSE)

plot(getverticeshr(d_GEPE_COPA_UD, percent = 95))

sqnm = d_GEPE_COPA_UD_area * 0.29155334959812 
sqnm


#put it through a speed filter
# --- Previous steps ---
# 1. Flatten your sf object, extract coordinates, drop geometry
df <- d_GEPE_COPA %>%
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  dplyr::select(Longitude, Latitude, datetime, Deployment)

# 2. Create a trip object and apply speed filter
trip_obj <- trip::trip(df, c("datetime", "Deployment"))
filtered_trip <- trip::speedfilter(trip_obj, max.speed = 6)
d_GEPE_COPA$filtered = filtered_trip
d_GEPE_COPA = d_GEPE_COPA %>%
  dplyr::filter(filtered == TRUE,
                Longitude < -48 & Longitude > -75, 
                Latitude < -60.5 & Latitude > -67 )

# --- New step: Convert filtered trip object back to sf ---
# 3. Convert to dataframe
d_GEPE_COPA <- as.data.frame(d_GEPE_COPA) %>%
  sf::st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326)

plot(st_geometry(d_GEPE_COPA))

sf::st_write(d_GEPE_COPA, dsn = "./Hinke_telemetry/GEPE_COPA_filtered.shp", delete_dsn = TRUE)

####Cape Shirreff----
d_GEPE_CS = d %>%
  dplyr::filter(Spp == "GEPE" & Site == "CS") %>% #REMOVE the single deployment responsible for the long trip into the Drake
  #dplyr::filter(month_day >= "11-15" | month_day <= "03-10") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(3031)

plot(st_geometry(d_GEPE_CS))

d_GEPE_CS_sp = as(d_GEPE_CS, "Spatial")

d_GEPE_CS_UD = adehabitatHR::kernelUD(d_GEPE_CS_sp, 
                                        h = "href", 
                                        grid = 500, 
                                        kern = "bivnorm", 
                                        same4all = TRUE)
d_GEPE_CS_UD_vert = getverticeshr(d_GEPE_CS_UD, percent = 95)
plot(d_GEPE_CS_UD_vert)

d_GEPE_CS_UD_area = kernel.area(d_GEPE_CS_UD, 
                                  percent = c(90),
                                  unin = c("m"),
                                  unout = c("km2"), standardize = FALSE)

plot(getverticeshr(d_GEPE_CS_UD, percent = 95))

sqnm = d_GEPE_CS_UD_area * 0.29155334959812 
sqnm

#ADELIES NOW
#Adelies
####Copacabana----
useless_months = c("12", "01")
d_ADPE_COPA = d %>%
  dplyr::filter(Spp == "ADPE" & Site == "COPA" & Month %in% useless_months & Longitude > -100 & Longitude < 90 & Deployment != "124") %>% 
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(3031)

plot(st_geometry(d_ADPE_COPA))

d_ADPE_COPA_sp = as(d_ADPE_COPA, "Spatial")

d_ADPE_COPA_UD = adehabitatHR::kernelUD(d_ADPE_COPA_sp, 
                                      h = "href", 
                                      grid = 500, 
                                      kern = "bivnorm", 
                                      same4all = TRUE)
d_ADPE_COPA_UD_vert = getverticeshr(d_ADPE_COPA_UD, percent = 95)
plot(d_ADPE_COPA_UD_vert)

d_ADPE_COPA_UD_area = kernel.area(d_ADPE_COPA_UD, 
                                percent = c(95),
                                unin = c("m"),
                                unout = c("km2"), standardize = FALSE)

plot(getverticeshr(d_ADPE_COPA_UD, percent = 95), add=TRUE)

sqnm = d_ADPE_COPA_UD_area * 0.29155334959812 
sqnm

# 2. Create a trip object and apply speed filter
trip_obj <- trip::trip(df, c("datetime", "Deployment"))
filtered_trip <- trip::speedfilter(trip_obj, max.speed = 6)
d_ADPE_COPA$filtered = filtered_trip
d_ADPE_COPA = d_ADPE_COPA %>%
  dplyr::filter(filtered == TRUE,
                Longitude < -48 & Longitude > -75, 
                Latitude < -60.5 & Latitude > -72 )

# --- New step: Convert filtered trip object back to sf ---
# 3. Convert to dataframe
d_ADPE_COPA <- as.data.frame(d_ADPE_COPA) %>%
  sf::st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326)

plot(st_geometry(d_ADPE_COPA))

sf::st_write(d_ADPE_COPA, dsn = "./Hinke_telemetry/ADPE_COPA_filtered.shp", delete_dsn = TRUE)



####FINALLY FUR SEALS----
d_AFS_CS = d %>%
  dplyr::filter(Spp == "AFS" & Site == "CS") %>%
  #dplyr::filter(month_day >= "11-15" | month_day <= "03-10") %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

plot(st_geometry(d_AFS_CS))

#put it through a speed filter
# --- Previous steps ---
# 1. Flatten your sf object, extract coordinates, drop geometry
df <- d_AFS_CS %>%
  mutate(x = sf::st_coordinates(.)[,1],
         y = sf::st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  dplyr::select(Longitude, Latitude, datetime, Deployment)

# 2. Create a trip object and apply speed filter
trip_obj <- trip::trip(df, c("datetime", "Deployment"))
filtered_trip <- trip::speedfilter(trip_obj, max.speed = 6)
d_AFS_CS$filtered = filtered_trip
d_AFS_CS = d_AFS_CS %>%
  dplyr::filter(filtered == TRUE,
                Longitude < -48, 
                Latitude < -50.5 & Latitude > -72 )

# --- New step: Convert filtered trip object back to sf ---
# 3. Convert to dataframe
d_AFS_CS <- as.data.frame(d_AFS_CS) %>%
  sf::st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326)

plot(st_geometry(d_AFS_CS))

sf::st_write(d_AFS_CS, dsn = "./Hinke_telemetry/AFS_CS_filtered.shp", delete_dsn = TRUE)