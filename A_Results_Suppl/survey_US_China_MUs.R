### Estimating LKB data for SSIW and BS MUs. 
## Using US and Chinese data [Wang et al. 2025. ICES J Mar Sci, 82(8)]
## Density extrapolated to MUs area

library(dplyr)

### --- Clean & aggregate US survey data ---
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)
# Filter low transect coverage (<10th percentil)
survey_filtered <- survey 
percentiles <- survey_filtered %>%
  group_by(gSSMU) %>%
  summarise(p10 = quantile(nmi.count, probs = 0.10, na.rm = TRUE))
percentiles$p10 <- ceiling(percentiles$p10/10)*10
survey_low_nmi <- survey_filtered %>%
  left_join(percentiles, by = "gSSMU") %>%
  filter(nmi.count >= p10)
survey <- survey_low_nmi
rm(percentiles, survey_filtered, survey_low_nmi)
### END changes
# # could try changing "biomass" in following line to "mean.density.gm2" or "median.density.gm2" but haven't done that
survey<-tapply(survey$mean.density.gm2,list(survey$Year,survey$gSSMU),mean,na.rm=TRUE)
survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                   gSSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                   density=c(survey),stringsAsFactors = FALSE)
survey$season<-ifelse(survey$cal.yr<2012,"S","W")
# use next line if want to remove winter survey data altogether (comment out if not desired)
#survey<-survey[survey$season=="S",]
survey$matchme<-paste(survey$cal.yr,survey$season,survey$gSSMU,sep="|")
# #print(str(survey))
survey <- survey %>% 
  filter(gSSMU %in% c("1", "2")) %>%
  mutate(
    cal.yr =  as.numeric(cal.yr),
    MU = ifelse(gSSMU == "1", "BS", "SSIW")
  ) %>%
  select(cal.yr, season, density, MU)


# Adding Chinese survey data
# Density data available for SSMUs APDPE, APDPW, APBSE, APBSW
# Fishing catches needs to be adjusted for same areas

# Create dataframe 'Chinese_north' (gSSMU2)
Chinese_SSIW <- data.frame(
  cal.yr  = c(2013, 2015, 2016, 2018, 2019),
  season  = rep("S", 5),
  density = c(18.0, 19.6, 38.5, 77.0, 62.9),  # g/m2
  MU   = rep("SSIW", 5),
  stringsAsFactors = FALSE)

Chinese_BS <- data.frame(
  cal.yr  = c(2013, 2015, 2016, 2018, 2019),
  season  = rep("S", 5),
  density = c(25.1, 17.6, 46.2, 35.6, 77.2),  # g/m2
  MU   = rep("BS", 5),
  stringsAsFactors = FALSE)

Chinese_survey <- rbind(Chinese_SSIW, Chinese_BS) 
rm(Chinese_SSIW, Chinese_BS)


### Combine and estimate MU biomass  ###
SSIW.km2 <- c(47134)
BS.km2 <- c(35208)

LKB <- rbind(survey, Chinese_survey)

LKB <- LKB %>% filter(!is.na(density)) %>%
  mutate(
    biomass = case_when(
      MU == "BS" ~ density*BS.km2,
      MU == "SSIW" ~ density*SSIW.km2))
# LKB are higher than in Wang et al. (2025), because they did not extrapolate to SSMU area.

### Export
openxlsx::write.xlsx(
  x = LKB,
  file = "./javier_analysis/Krill_data/LKB_MUs.xlsx",
  asTable = TRUE
)
