# Imputing LKB data based on SAM

library(dplyr)

# Import SAM data
sam<-read.csv("./Supplementary Files/sam.csv")
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$cal.yr<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$cal.yr,sam$season),mean)
sam<-data.frame(cal.yr=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))

# Formating
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
LKB_imputed <- rbind(LKB_imputed_gSSM1, LKB_imputed_gSSM2)
rm(LKB_imputed_gSSM1,LKB_imputed_gSSM2)
