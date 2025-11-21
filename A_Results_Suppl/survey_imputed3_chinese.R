### Imputing LKB data based on SAM + Wang et al. 2025. ICES J Mar Sci, 82(8)
## Density extrapolated to SSMUs area

library(dplyr)

# Import SAM data
sam<-read.csv("./Supplementary Files/sam.csv")
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$cal.yr<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$cal.yr,sam$season),mean)
sam<-data.frame(cal.yr=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))

# Formatting
sam <- sam %>% mutate(cal.yr = as.integer(cal.yr))
sam <- sam %>% filter(season == "S") %>%   # only imputed summer data
  filter(cal.yr > 1979) %>%                # remove data before 1980, and between 1996-2011 (when surveys occurred)
  filter(cal.yr < 1996 | cal.yr > 2011)

LKB_imputed_gSSM1 <- sam %>%
  mutate(LKB = ifelse(sam<0, 1291502, 1818980),
         gSSMU = "1")
LKB_imputed_gSSM2 <- sam %>%
  mutate(LKB = ifelse(sam<0, 3045776, 7265616),
         gSSMU = "2")
LKB_imputed2 <- rbind(LKB_imputed_gSSM1, LKB_imputed_gSSM2)
rm(LKB_imputed_gSSM1,LKB_imputed_gSSM2)

# Adding Chinese survey data
# Density data available for SSMUs APDPE, APDPW, APBSE, APBSW
# Fishing catches needs to be adjusted for same areas

# Create dataframe 'Chinese_north' (gSSMU2)
Chinese_SSIW <- data.frame(
  cal.yr  = c(2013, 2015, 2016, 2018, 2019),
  season  = rep("S", 5),
  LKB1    = c(418000, 400000, 896000, 1794000, 1186000),
  density = c(18.0, 19.6, 38.5, 77.0, 62.9),  # g/m2
  gSSMU   = rep("2", 5),
  stringsAsFactors = FALSE)

Chinese_BS <- data.frame(
  cal.yr  = c(2013, 2015, 2016, 2018, 2019),
  season  = rep("S", 5),
  LKB1    = c(563000, 395000, 1035000, 798000, 1730000),
  density = c(25.1, 17.6, 46.2, 35.6, 77.2),  # g/m2
  gSSMU   = rep("1", 5),
  stringsAsFactors = FALSE)

Chinese_survey <- rbind(Chinese_SSIW, Chinese_BS) 
rm(Chinese_SSIW, Chinese_BS)
### Run only first time
# # Obtain area for each SSMU:
# library(CCAMLRGIS)
# SSMU = load_SSMUs()
# SSMU <- SSMU %>% filter(GAR_Short_Label %in% c("APDPE", "APDPW", "APBSE", "APBSW"))
# ## Clip out land
# #Load Coastline
# Coast=load_Coastline()
# 
# #Isolate land and merge (union) polygons into one:
# Land=Coast[Coast$surface=="Land",]
# Land=st_union(Land)
# 
# #Clip polygons
# SSMU_pelagic = suppressWarnings(st_difference(SSMU,Land))
# 
# #Estimate Surface Area (Km2)
# Ar=round(st_area(SSMU_pelagic)/1000000,1)
# SSMU_pelagic$AreaKm2=as.numeric(Ar)
# 
# # Summarize
# SSMU_area <- SSMU_pelagic %>% select(GAR_Short_Label, AreaKm2) 

### SSMU area are:
SSMUkm2 <- data.frame(
  SSMU = c("APDPW", "APDPE", "APBSW", "APBSE"),
  Area_Km2  = c(15956, 16554, 22155, 28815))

### Estimate Biomass for each gSSMU
area1 <- SSMUkm2$Area_Km2[3] + SSMUkm2$Area_Km2[4]
area2 <- SSMUkm2$Area_Km2[1] + SSMUkm2$Area_Km2[2]

Chinese_survey <- Chinese_survey %>%
  mutate(
    LKB2 = case_when(
      gSSMU == "1" ~ density*area1,
      gSSMU == "2" ~ density*area2))
# LKB are higher than in Wang et al. (2025), because they did not extrapolate to SSMU area.

### Add Chinese survey data into LKB imputed dataframe
LKB_imputed_Chinese <- left_join(LKB_imputed2, Chinese_survey,
                    by = c("cal.yr", "gSSMU"))

## Chinese 1: original data LKB1 -> LKB
LKB_imputed_Chinese$LKB[!is.na(LKB_imputed_Chinese$LKB1)] <- 
  LKB_imputed_Chinese$LKB1[!is.na(LKB_imputed_Chinese$LKB1)]
LKB_imputed_Chinese1 <- LKB_imputed_Chinese %>%
  select(cal.yr, gSSMU, season.x, LKB) %>%
  rename(season = season.x)

## Chinese 2: extrapolated data LKB2 -> LKB
LKB_imputed_Chinese$LKB[!is.na(LKB_imputed_Chinese$LKB2)] <- 
  LKB_imputed_Chinese$LKB2[!is.na(LKB_imputed_Chinese$LKB2)]
LKB_imputed_Chinese2 <- LKB_imputed_Chinese %>%
  select(cal.yr, gSSMU, season.x, LKB) %>%
  rename(season = season.x)

### end

