### Watter et al. 2020
# Model imputed LKB

library(dplyr)
library(tidyverse)


#### First, import data from the model ####

# krill survey biomass
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)
survey<-tapply(survey$biomass,list(survey$Year,survey$gSSMU),mean,na.rm=TRUE)
survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                   gSSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                   survey=c(survey),stringsAsFactors = FALSE)
survey$season<-ifelse(survey$cal.yr<2012,"S","W")
survey$matchme<-paste(survey$cal.yr,survey$season,survey$gSSMU,sep="|")

survey <- survey %>% mutate(cal.yr = as.integer(cal.yr))
survey <- survey %>% filter(season == "S")   # only imputed summer data
survey <- survey %>% mutate(gSSMU = as.integer(gSSMU))

# Import SAM data
sam<-read.csv("./Supplementary Files/sam.csv")
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$cal.yr<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$cal.yr,sam$season),mean)
sam<-data.frame(cal.yr=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))

sam <- sam %>% mutate(cal.yr = as.integer(cal.yr))
sam <- sam %>% filter(season == "S")   # only imputed summer data
  

# krill fishery catches
fishery<-read.csv("./Supplementary Files/c1.csv",header=TRUE,stringsAsFactors = FALSE)
fishery$season<-ifelse(is.element(fishery$Month,c(10:12,1:3)),"S","W")
gSSMU1<-c("APBSE","APBSW")
gSSMU2<-c("APDPE","APDPW","APEI")
gSSMU3<-"APPA"
gSSMU4<-c("APW","APE")
fishery$gSSMU<-rep(NA,dim(fishery)[1])
fishery$gSSMU<-ifelse(is.element(fishery$AssignedSSMU,gSSMU1),1,fishery$gSSMU)
fishery$gSSMU<-ifelse(is.element(fishery$AssignedSSMU,gSSMU2),2,fishery$gSSMU)
fishery$gSSMU<-ifelse(is.element(fishery$AssignedSSMU,gSSMU3),3,fishery$gSSMU)
fishery$gSSMU<-ifelse(is.element(fishery$AssignedSSMU,gSSMU4),4,fishery$gSSMU)
fishery<-fishery[!is.na(fishery$gSSMU),]
fishery<-tapply(fishery$TotalCatch,list(fishery$CalendarYear,fishery$gSSMU,fishery$season),sum)
fishery<-data.frame(cal.yr=rep(dimnames(fishery)[[1]],dim(fishery)[2]*dim(fishery)[3]),
                    gSSMU=rep(rep(dimnames(fishery)[[2]],each=dim(fishery)[1]),dim(fishery)[3]),
                    season=rep(dimnames(fishery)[[3]],each=dim(fishery)[1]*dim(fishery)[2]),
                    catch=c(fishery),stringsAsFactors = FALSE)
fishery$cal.yr<-as.numeric(as.character(fishery$cal.yr))
fishery$gSSMU<-as.numeric(as.character(fishery$gSSMU))
fishery$matchme<-paste(fishery$cal.yr,fishery$season,fishery$gSSMU,sep="|")

fishery <- fishery %>% mutate(cal.yr = as.integer(cal.yr)) 
fishery <- fishery %>% filter(season == "S")   # only imputed summer data

# junk <- readRDS("./javier_analysis/Watters_model_output1/junk.rds")
# data <- junk %>% filter(season == "S") %>% filter(!is.na(survey)) %>% distinct(survey)

#### Imputing LKB data based on SAM ####

# join variables
data <- left_join(fishery, survey, by = c("cal.yr", "gSSMU", "season", "matchme")) 
data <- data %>% filter(gSSMU %in% c(1,2))
data <- left_join(data, sam, by = c("cal.yr", "season")) 
data <- data %>% mutate(sam.sign = (ifelse(sam<0, "Neg", "Pos")))

# Estimating mean(logLKB) for survey data
data$LnSurvey <- log(data$survey)

meanLogS <- data.frame(SAM.pos = numeric(2), SAM.neg = numeric(2))
sdLogS <- data.frame(SAM.pos = numeric(2), SAM.neg = numeric(2))

data.pos <- data %>% filter(sam.sign == "Pos")
data.neg <- data %>% filter(sam.sign == "Neg")

for(i in 1:2) {
  gS <- data.neg %>% filter(!is.na(survey)) %>% 
    filter(gSSMU == i) %>% filter(season == "S")
  meanLogS[i,1] <- mean(gS$LnSurvey)
  sdLogS[i,1] <- sd(gS$LnSurvey)
}
for(i in 1:2) {
  gS <- data.pos %>% filter(!is.na(survey)) %>% 
    filter(gSSMU == i) %>% filter(season == "S")
  meanLogS[i,2] <- mean(gS$LnSurvey)
  sdLogS[i,2] <- sd(gS$LnSurvey)
}


# Imputing missing data
data <- data %>%
  mutate(
    # 0 if survey available, 1 if missing
    impute.me = if_else(is.na(survey), 1, 0),
    
    # fill log survey
    LogSurvey_imputed = case_when(
      impute.me == 1L & gSSMU == 1 ~ meanLogS[1, "mean"],
      impute.me == 1L & gSSMU == 2 ~ meanLogS[2, "mean"],
      impute.me == 0L              ~ LnSurvey,          # or log(survey)
      TRUE                         ~ NA_real_
    ),
    
    # back-transform
    survey_imputed = exp(LogSurvey_imputed)
  )


mulogsummer <- matrix(NA_real_, nrow = 2, ncol = 2)

mulogsummer

K <- data.frame()

## ---- Priors -----------------------------------------------------------
# mulogsummer[i,j] ~ dunif(0.1 * meanlogsummer[i,j], 10 * meanlogsummer[i,j])
mulogsummer <- matrix(NA_real_, nrow = 2, ncol = 2)
for (i in 1:2) {
  for (j in 1:2) {
    mulogsummer[i, j] <- runif(
      n   = 1,
      min = 0.1 * meanlogsummer[i, j],
      max = 10  * meanlogsummer[i, j]
    )
  }
}

# sigmalogsummer ~ dunif(0.1 * sdlogsummer, 10 * sdlogsummer)
sigmalogsummer <- runif(
  n   = 1,
  min = 0.1 * sdlogsummer,
  max = 10  * sdlogsummer
)
taulogsummer <- sigmalogsummer^(-2)



for(i in nrow(missLKB)) {
  
}


### Waters original Model script ####

# likelihood
for(i in 1:nsummerobs){
  # I think this structure will only work if the data are arranged such that all summer data are first followed by winter data
  # have never observed gSSMU-scale biomass < 10,000t or > 100,000,000t
  # catches at gSSMU scale have ranged from 0t to about 117,000t (during period of study)
  # assume catch has never been greater than biomass so truncate distribution of summer biomasses as follows
  # lower bound of max(10000t, catch in year i) and upper bound of 100000000t
  lower[i]<-max(10000,catch[i])
  summer[i]~dlnorm(mulogsummer[gssmu[i],samclass[i]],taulogsummer) T(lower[i],100000000)
}

# priors
# this prior specification patterned after http://doingbayesiandataanalysis.blogspot.com/2016/04/bayesian-estimation-of-log-normal.html
for(i in 1:2){ # two gSSMUs
  for(j in 1:2){ # two SAM classes
    mulogsummer[i,j]~dunif(0.1*meanlogsummer[i,j],10*meanlogsummer[i,j])
  }
}
taulogsummer<-pow(sigmalogsummer,-2)
sigmalogsummer~dunif(0.1*sdlogsummer,10*sdlogsummer)

# substitute imputed harvest rates
for(i in 1:nsummerobs){
  hr.summer[i]<-ifelse(impute.me[i]==1,catch[i]/summer[i],1)
  bmass.summer[i]<-ifelse(impute.me[i]==1,summer[i],1)
}
for(i in (nsummerobs+1):nobs){
  hr.summer[i]<-0
  bmass.summer[i]<-0
}

for(i in 1:nobs){
  hr[i]<-ifelse(impute.me[i]==1,hr.summer[i],catch[i]/survey[i])
  hrclass[i]<-ifelse(hr[i]<=0.01,1,ifelse(hr[i]>=0.1,3,2))
  bmass[i]<-ifelse(impute.me[i]==1,bmass.summer[i],survey[i])
  bclass[i]<-ifelse(bmass[i]<=1000000,1,2)
}
