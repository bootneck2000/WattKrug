
library(sf)
library(terra)
library(tidyverse)


whales_1 = terra::rast("/Users/god/Documents/R workspace/WattKrug/zenodo_Whales/Cruise 1.grd")
whales_2 = terra::rast("/Users/god/Documents/R workspace/WattKrug/zenodo_Whales/Cruise 2.grd")
whales_3 = terra::rast("/Users/god/Documents/R workspace/WattKrug/zenodo_Whales/Cruise 3.grd")

whales_1 <- ifel(is.na(whales_1), 0, whales_1)
whales_2 <- ifel(is.na(whales_2), 0, whales_2)
whales_3 <- ifel(is.na(whales_3), 0, whales_3)

ssmu = sf::read_sf("/Volumes/GodDrive/Quantarctica3/EnvironmentalManagement/CCAMLR/CCAMLR_SSMUs.shp") %>%
  sf::st_transform(crs = sf::st_crs(whales_1))
####create objects for each SSMU, and for each gSSMU----
bsw = ssmu[5, ] %>%
  terra::vect() 
bse = ssmu[6,] %>%
  terra::vect()
dpe = ssmu[4, ] %>%
  terra::vect() 
dpw = ssmu[3,] %>%
  terra::vect()
apei = ssmu[7, ] %>%
  terra::vect() 

gSSMU_1 = rbind(bsw, bse) %>%
  terra::aggregate()

gSSMU_2 = rbind(dpw,dpe,apei) %>%
  terra::aggregate()

#humpback whale consumptions by SSMU
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

#Now overlap between gSSMU in % by penguin UD - set CRS first
whale_proj = terra::crs(whales_1, proj = TRUE)
####DPW----
dpw_whales_1 = terra::extract(whales_1, dpw, fun = sum)
dpw_whales_2 = terra::extract(whales_2, dpw, fun = sum)
dpw_whales_3 = terra::extract(whales_3, dpw, fun = sum)

dpw_whales_1
dpw_whales_2
dpw_whales_3

dpw_whales_1_consumption_day = dpw_whales_1$z * 0.632
dpw_whales_2_consumption_day = dpw_whales_2$z * 0.632
dpw_whales_3_consumption_day = dpw_whales_3$z * 0.632

dpw_whales_1_breed_season = dpw_whales_1_consumption_day * 30
dpw_whales_2_breed_season = dpw_whales_2_consumption_day * 30
dpw_whales_3_breed_season = dpw_whales_3_consumption_day * 30

dpw_whales_total_breed_season = dpw_whales_1_breed_season + dpw_whales_2_breed_season + dpw_whales_3_breed_season

####DPE----
dpe_whales_1 = terra::extract(whales_1, dpe, fun = sum)
dpe_whales_2 = terra::extract(whales_2, dpe, fun = sum)
dpe_whales_3 = terra::extract(whales_3, dpe, fun = sum)

dpe_whales_1
dpe_whales_2
dpe_whales_3

dpe_whales_1_consumption_day = dpe_whales_1$z * 0.632
dpe_whales_2_consumption_day = dpe_whales_2$z * 0.632
dpe_whales_3_consumption_day = dpe_whales_3$z * 0.632

dpe_whales_1_breed_season = dpe_whales_1_consumption_day * 30
dpe_whales_2_breed_season = dpe_whales_2_consumption_day * 30
dpe_whales_3_breed_season = dpe_whales_3_consumption_day * 30

dpe_whales_total_breed_season = dpe_whales_1_breed_season + dpe_whales_2_breed_season + dpe_whales_3_breed_season

####APEI----
apei_whales_1 = terra::extract(whales_1, apei, fun = sum)
apei_whales_2 = terra::extract(whales_2, apei, fun = sum)
apei_whales_3 = terra::extract(whales_3, apei, fun = sum)

apei_whales_1
apei_whales_2
apei_whales_3

apei_whales_1_consumption_day = apei_whales_1$z * 0.632
apei_whales_2_consumption_day = apei_whales_2$z * 0.632
apei_whales_3_consumption_day = apei_whales_3$z * 0.632

apei_whales_1_breed_season = apei_whales_1_consumption_day * 30
apei_whales_2_breed_season = apei_whales_2_consumption_day * 30
apei_whales_3_breed_season = apei_whales_3_consumption_day * 30

apei_whales_total_breed_season = apei_whales_1_breed_season + apei_whales_2_breed_season + apei_whales_3_breed_season

#Penguins now

#humpback whale consumption by penguin UD
#####Cape Shireff Chinstraps----
d_CHPE_CS_trans_UD = d_CHPE_CS %>%
  sf::st_transform(crs = whale_proj) %>%
  sf::as_Spatial() %>%
  adehabitatHR::kernelUD(h = "href", 
                         grid = 500, 
                         kern = "bivnorm", 
                         same4all = TRUE) %>%
  adehabitatHR::getverticeshr(percent = 95) %>%
  terra::vect()
#%overlap between pengo MCP and SSMU / gSSMU
ov <- terra::intersect(d_CHPE_CS_trans_UD, gSSMU_2)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_2,  unit = "m"), na.rm = TRUE)  # total area of y

d_CHPE_CS_pct_y_within_x <- 100 * (a_xy / a_y)
d_CHPE_CS_pct_y_within_x

#consumption within the UD
bsw_whales_1_d_CHPE_CS = terra::extract(whales_1, d_CHPE_CS_trans_UD, fun = sum)
bsw_whales_2_d_CHPE_CS = terra::extract(whales_2, d_CHPE_CS_trans_UD, fun = sum)
bsw_whales_3_d_CHPE_CS = terra::extract(whales_3, d_CHPE_CS_trans_UD, fun = sum)

bsw_whales_1_d_CHPE_CS_consumption_day = bsw_whales_1_d_CHPE_CS$z * 0.632
bsw_whales_2_d_CHPE_CS_consumption_day = bsw_whales_2_d_CHPE_CS$z * 0.632
bsw_whales_3_d_CHPE_CS_consumption_day = bsw_whales_3_d_CHPE_CS$z * 0.632

bsw_whales_1_breed_season_d_CHPE_CS = bsw_whales_1_d_CHPE_CS_consumption_day * 30
bsw_whales_2_breed_season_d_CHPE_CS = bsw_whales_2_d_CHPE_CS_consumption_day * 30
bsw_whales_3_breed_season_d_CHPE_CS = bsw_whales_3_d_CHPE_CS_consumption_day * 30

bsw_whales_total_breed_season_d_CHPE_CS = bsw_whales_1_breed_season_d_CHPE_CS + bsw_whales_2_breed_season_d_CHPE_CS + bsw_whales_3_breed_season_d_CHPE_CS


####Copacabana Chinstraps----
d_CHPE_COPA_trans_UD = d_CHPE_COPA %>%
  sf::st_transform(crs = whale_proj) %>%
  sf::as_Spatial() %>%
  adehabitatHR::kernelUD(h = "href", 
                         grid = 500, 
                         kern = "bivnorm", 
                         same4all = TRUE) %>%
  adehabitatHR::getverticeshr(percent = 95) %>%
  terra::vect()
#%overlap between pengo UD and SSMU / gSSMU
ov <- terra::intersect(d_CHPE_COPA_trans_UD, gSSMU_1)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_1,  unit = "m"), na.rm = TRUE)  # total area of y

d_CHPE_COPA_pct_y_within_x <- 100 * (a_xy / a_y)
d_CHPE_COPA_pct_y_within_x

#consumption within the UD
bsw_whales_1_d_CHPE_COPA = terra::extract(whales_1, d_CHPE_COPA_trans_UD, fun = sum)
bsw_whales_2_d_CHPE_COPA = terra::extract(whales_2, d_CHPE_COPA_trans_UD, fun = sum)
bsw_whales_3_d_CHPE_COPA = terra::extract(whales_3, d_CHPE_COPA_trans_UD, fun = sum)

bsw_whales_1_d_CHPE_COPA_consumption_day = bsw_whales_1_d_CHPE_COPA$z * 0.632
bsw_whales_2_d_CHPE_COPA_consumption_day = bsw_whales_2_d_CHPE_COPA$z * 0.632
bsw_whales_3_d_CHPE_COPA_consumption_day = bsw_whales_3_d_CHPE_COPA$z * 0.632

bsw_whales_1_breed_season_d_CHPE_COPA = bsw_whales_1_d_CHPE_COPA_consumption_day * 30
bsw_whales_2_breed_season_d_CHPE_COPA = bsw_whales_2_d_CHPE_COPA_consumption_day * 30
bsw_whales_3_breed_season_d_CHPE_COPA = bsw_whales_3_d_CHPE_COPA_consumption_day * 30

bsw_whales_total_breed_season_d_CHPE_COPA = bsw_whales_1_breed_season_d_CHPE_COPA + bsw_whales_2_breed_season_d_CHPE_COPA + bsw_whales_3_breed_season_d_CHPE_COPA

####Copacabana Adelies

####Copacabana Adelies----
d_ADPE_COPA_trans_UD = d_ADPE_COPA %>%
  sf::st_transform(crs = whale_proj) %>%
  sf::as_Spatial() %>%
  adehabitatHR::kernelUD(h = "href", 
                         grid = 500, 
                         kern = "bivnorm", 
                         same4all = TRUE) %>%
  adehabitatHR::getverticeshr(percent = 95) %>%
  terra::vect()
#%overlap between pengo UD and SSMU / gSSMU
ov <- terra::intersect(d_ADPE_COPA_trans_UD, gSSMU_1)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_1,  unit = "m"), na.rm = TRUE)  # total area of y

d_ADPE_COPA_pct_y_within_x <- 100 * (a_xy / a_y)
d_ADPE_COPA_pct_y_within_x

#consumption within the UD
bsw_whales_1_d_ADPE_COPA = terra::extract(whales_1, d_ADPE_COPA_trans_UD, fun = sum)
bsw_whales_2_d_ADPE_COPA = terra::extract(whales_2, d_ADPE_COPA_trans_UD, fun = sum)
bsw_whales_3_d_ADPE_COPA = terra::extract(whales_3, d_ADPE_COPA_trans_UD, fun = sum)

bsw_whales_1_d_ADPE_COPA_consumption_day = bsw_whales_1_d_ADPE_COPA$z * 0.632
bsw_whales_2_d_ADPE_COPA_consumption_day = bsw_whales_2_d_ADPE_COPA$z * 0.632
bsw_whales_3_d_ADPE_COPA_consumption_day = bsw_whales_3_d_ADPE_COPA$z * 0.632

bsw_whales_1_breed_season_d_ADPE_COPA = bsw_whales_1_d_ADPE_COPA_consumption_day * 30
bsw_whales_2_breed_season_d_ADPE_COPA = bsw_whales_2_d_ADPE_COPA_consumption_day * 30
bsw_whales_3_breed_season_d_ADPE_COPA = bsw_whales_3_d_ADPE_COPA_consumption_day * 30

bsw_whales_total_breed_season_d_ADPE_COPA = bsw_whales_1_breed_season_d_ADPE_COPA + bsw_whales_2_breed_season_d_ADPE_COPA + bsw_whales_3_breed_season_d_ADPE_COPA

####Cape Shirreff Gentoos----
d_GEPE_CS_trans_UD = d_GEPE_CS %>%
  sf::st_transform(crs = whale_proj) %>%
  sf::as_Spatial() %>%
  adehabitatHR::kernelUD(h = "href", 
                         grid = 500, 
                         kern = "bivnorm", 
                         same4all = TRUE) %>%
  adehabitatHR::getverticeshr(percent = 95) %>%
  terra::vect()
#%overlap between pengo UD and SSMU / gSSMU
ov <- terra::intersect(d_GEPE_CS_trans_UD, gSSMU_2)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_2,  unit = "m"), na.rm = TRUE)  # total area of y

d_GEPE_CS_pct_y_within_x <- 100 * (a_xy / a_y)
d_GEPE_CS_pct_y_within_x

#CS Gentoos also overlap with gSSMU #1, presumably during winter.  So we need to look at % overlap here too and add to table
ov <- terra::intersect(d_GEPE_CS_trans_UD, gSSMU_1)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_1,  unit = "m"), na.rm = TRUE)  # total area of y

d_GEPE_CS_pct_y_within_x_gssmu1 <- 100 * (a_xy / a_y)
d_GEPE_CS_pct_y_within_x_gssmu1

#consumption within the UD
bsw_whales_1_d_GEPE_CS = terra::extract(whales_1, d_GEPE_CS_trans_UD, fun = sum)
bsw_whales_2_d_GEPE_CS = terra::extract(whales_2, d_GEPE_CS_trans_UD, fun = sum)
bsw_whales_3_d_GEPE_CS = terra::extract(whales_3, d_GEPE_CS_trans_UD, fun = sum)

bsw_whales_1_d_GEPE_CS_consumption_day = bsw_whales_1_d_GEPE_CS$z * 0.632
bsw_whales_2_d_GEPE_CS_consumption_day = bsw_whales_2_d_GEPE_CS$z * 0.632
bsw_whales_3_d_GEPE_CS_consumption_day = bsw_whales_3_d_GEPE_CS$z * 0.632

bsw_whales_1_breed_season_d_GEPE_CS = bsw_whales_1_d_GEPE_CS_consumption_day * 30
bsw_whales_2_breed_season_d_GEPE_CS = bsw_whales_2_d_GEPE_CS_consumption_day * 30
bsw_whales_3_breed_season_d_GEPE_CS = bsw_whales_3_d_GEPE_CS_consumption_day * 30

bsw_whales_total_breed_season_d_GEPE_CS = bsw_whales_1_breed_season_d_GEPE_CS + bsw_whales_2_breed_season_d_GEPE_CS + bsw_whales_3_breed_season_d_GEPE_CS

####Copacabana Gentoos----
d_GEPE_COPA_trans_UD = d_GEPE_COPA %>%
  sf::st_transform(crs = whale_proj) %>%
  sf::as_Spatial() %>%
  adehabitatHR::kernelUD(h = "href", 
                         grid = 500, 
                         kern = "bivnorm", 
                         same4all = TRUE) %>%
  adehabitatHR::getverticeshr(percent = 95) %>%
  terra::vect()
#%overlap between pengo UD and SSMU / gSSMU
ov <- terra::intersect(d_GEPE_COPA_trans_UD, gSSMU_1)
a_xy <- sum(expanse(ov, unit = "m"), na.rm = TRUE)  # overlap area
a_y  <- sum(expanse(gSSMU_1,  unit = "m"), na.rm = TRUE)  # total area of y

d_GEPE_COPA_pct_y_within_x <- 100 * (a_xy / a_y)
d_GEPE_COPA_pct_y_within_x

#consumption within the UD
bsw_whales_1_d_GEPE_COPA = terra::extract(whales_1, d_GEPE_COPA_trans_UD, fun = sum)
bsw_whales_2_d_GEPE_COPA = terra::extract(whales_2, d_GEPE_COPA_trans_UD, fun = sum)
bsw_whales_3_d_GEPE_COPA = terra::extract(whales_3, d_GEPE_COPA_trans_UD, fun = sum)

bsw_whales_1_d_GEPE_COPA_consumption_day = bsw_whales_1_d_GEPE_COPA$z * 0.632
bsw_whales_2_d_GEPE_COPA_consumption_day = bsw_whales_2_d_GEPE_COPA$z * 0.632
bsw_whales_3_d_GEPE_COPA_consumption_day = bsw_whales_3_d_GEPE_COPA$z * 0.632

bsw_whales_1_breed_season_d_GEPE_COPA = bsw_whales_1_d_GEPE_COPA_consumption_day * 30
bsw_whales_2_breed_season_d_GEPE_COPA = bsw_whales_2_d_GEPE_COPA_consumption_day * 30
bsw_whales_3_breed_season_d_GEPE_COPA = bsw_whales_3_d_GEPE_COPA_consumption_day * 30

bsw_whales_total_breed_season_d_GEPE_COPA = bsw_whales_1_breed_season_d_GEPE_COPA + bsw_whales_2_breed_season_d_GEPE_COPA + bsw_whales_3_breed_season_d_GEPE_COPA

####Table for %overlap between penguin UDs and gSSMU/SSMU----
tab <- tibble::tibble(
  row = c("Adelie", "Chinstrap", "Gentoo"),
  
  # COPA (gSSMU #1)
  gSSMU_1 = c(d_ADPE_COPA_pct_y_within_x/100, d_CHPE_COPA_pct_y_within_x/100, d_GEPE_COPA_pct_y_within_x/100),
  gSSMU_2 = c(NA,  d_CHPE_CS_pct_y_within_x/100, d_GEPE_CS_pct_y_within_x/100),
)

tab_mat <- tab %>%
  tibble::column_to_rownames("row") %>%
  as.matrix()

tab_mat <- apply(tab_mat, 2, function(x) sprintf("%.2f", as.numeric(x)))
rownames(tab_mat) <- c("Adelie", "Chinstrap", "Gentoo")

extra_val = as.character(round((d_GEPE_CS_pct_y_within_x_gssmu1/100), digits = 2))
tab_mat[3, 1] <- paste0(tab_mat[3, 1], " (", extra_val, ")")

kbl_obj <-
  knitr::kable(
    tab_mat,
    format = "html",
    booktabs = TRUE,
    align = "c",
    caption = "Table SX. Overlap (%) of pygoscelid penguin 95% UD with the gSSMU they were assigned to.  Number in parentheses 
    indicates the % overlap of post-breeding tracking data for Gentoo penguins at Cape Shirreff, assigned to gSSMU_2 during summer but 
    move to gSSMU_1 after breeding, exclusively in SSMU Bransfield Strait West."
  ) %>%
  kableExtra::kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "condensed"),
    position = "center",
    font_size = 12
  ) %>%
  # (Optional) emphasize row names column
  kableExtra::column_spec(1, bold = TRUE, width = "1.5in") %>%
  # (Optional) set widths for value columns
  kableExtra::column_spec(2:3, width = "1.5in") %>%
  kableExtra::column_spec(1, extra_css = "text-align: right;")
kbl_obj
#Fisheries removals by SSMU
####Now we read in the catch data from Watters et al, subset the two SSMU, sum the catch for the same periods, and compare----
summer_months = c(11, 12, 1, 2)
winter_months = c(3, 4, 5, 6, 7, 8, 9, 10)
catch_BSW = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APBSW") & CalendarYear >= "1982" & CalendarYear <= "2016" & Month %in% summer_months) %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

total_catch_BSW = sum (catch_BSW$TotalCatch, na.rm = TRUE)

catch_BSE = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APBSE") & CalendarYear >= "1982" & CalendarYear <= "2016" & Month %in% summer_months)%>%#, CalendarYear >= "2010") %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

total_catch_BSE = sum (catch_BSW$TotalCatch, na.rm = TRUE)

catch_DPW = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APDPW") & CalendarYear >= "1982" & CalendarYear <= "2016" & Month %in% summer_months) %>%#, CalendarYear >= "2010") %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

total_catch_DPW = sum (catch_BSW$TotalCatch, na.rm = TRUE)

catch_DPE = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APDPE") & CalendarYear >= "1982" & CalendarYear <= "2016" & Month %in% summer_months)%>%#, CalendarYear >= "2010") %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

total_catch_DPE = sum (catch_BSW$TotalCatch, na.rm = TRUE)

catch_EI = read_csv("/Users/god/Documents/R workspace/WattKrug/Supplementary Files/c1.csv") %>%
  dplyr::filter(AssignedSSMU == c("APEI") & CalendarYear >= "1982" & CalendarYear <= "2016" & Month %in% summer_months) %>%#, CalendarYear >= "2010") %>%
  dplyr::group_by(AssignedSSMU, CalendarYear) %>%
  dplyr::summarise(TotalCatch = sum(TotalCatch, na.rm = TRUE))

total_catch_EI = sum (catch_BSW$TotalCatch, na.rm = TRUE)

####total catches by gSSMU for whales and fishery----
total_catch_gssmu_1 = total_catch_BSW + total_catch_BSE
total_catch_gssmu_2 = total_catch_DPW + total_catch_DPE + total_catch_EI

gssmu1_whales_1 = terra::extract(whales_1, gSSMU_1, fun = sum)
gssmu1_whales_2 = terra::extract(whales_2, gSSMU_1, fun = sum)
gssmu1_whales_3 = terra::extract(whales_3, gSSMU_1, fun = sum)

x = c(gssmu1_whales_1$z,gssmu1_whales_2$z,gssmu1_whales_3$z)
(mean(x) *0.632) *90
sd(x)

####back of the fag packt reverse calc for male AFS to hit the 155kt 48.1 trigger----
#uses Boyd et al. consumption estimate of a 4yr old male eating 3tonnes/year, or an average of 8kg/day

quota = 31000
seal_consumption = 0.008*90 #8kg per day for 90days, or between January 1st and March 31st each year
number_seals = (quota/seal_consumption)
number_seals


catch_and_seals = catch %>%
  dplyr::mutate(equivalent_AFS_n = TotalCatch/0.72)