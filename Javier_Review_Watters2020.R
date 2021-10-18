# Analisis Watters et al 2020


library(readr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(tidyverse)
library(stringr)


S1 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM2_ESM.csv") # Krill catches
S2 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM3_ESM.csv") # Clutch dates
S3 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM4_ESM.csv", col_types="ccccdddd") # Egg mass and volume
S4 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM5_ESM.csv") # Fledgling mass
S5 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM6_ESM.csv") # LKB
S6 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM7_ESM.csv") # Adult penguin mass at laying
S7 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM8_ESM.csv") # ONI
S8 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM9_ESM.csv") # Cohort strength
S9 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM10_ESM.csv") # SAM
S10 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM11_ESM.csv") # Chick breeding success
S11 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM12_ESM.csv") # Forging trip duration
S12 <- read_csv("./Supplementary Files/41598_2020_59223_MOESM13_ESM.csv") # Krill catches and LHR



###############################
# Relationship SAM and LKB ------------------------------------------------
###############################

# LKB vs ONI & SAM
#Estimate Summer and Winter Indexes - ONI
oni <- S7
oni$yr<-ifelse(is.element(oni$SEAS,c("OND","NDJ")),oni$YR+1,oni$YR)
oni$season<-ifelse(is.element(oni$SEAS,c("OND","NDJ","DJF","JFM")),"S",NA)
oni$season<-ifelse(is.element(oni$SEAS,c("AMJ","MJJ","JJA","JAS")),"W",oni$season)
oni<-na.omit(oni)
oni<-tapply(oni$ANOM,list(oni$yr,oni$season),mean)
oni<-data.frame(yr=rep(dimnames(oni)[[1]],2),season=rep(dimnames(oni)[[2]],each=dim(oni)[1]),oni=c(oni))
names(oni)<-c("Year","season","oni")


#Estimate Summer and Winter Indexes - SAM
sam <- S9
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$Year<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$Year,sam$season),mean)
sam<-data.frame(Year=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))


# LKB 
survey<-S5[S5$nmi.count>=50,] # remove extreme bad surveys
survey$biomass <- as.integer(survey$biomass)
survey<-tapply(survey$biomass,list(survey$Year,survey$gSSMU),mean,na.rm=TRUE)
survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                   gSSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                   survey=c(survey),stringsAsFactors = FALSE)
survey$season<-ifelse(survey$cal.yr<2012,"S","W")
survey<-na.omit(survey)
names(survey) <- c("Year","gSSMU","biomass", "season")


# LKB vs SAM
LKB <- merge(x = survey, y = sam, by = "Year", all.x = TRUE)
LKB <- merge(x = LKB, y = oni, by = "Year", all.x = TRUE)
names(LKB)    
names(LKB) <- c("Year",  "gSSMU",  "biomass",  "season.b", "season.sam", "sam",  "season.oni", "oni")
LKB$Year <- as.integer(LKB$Year)

# Watters et al. 2020: "We modeled LKB as a function of stratum and the sign of the SAM during summer." 
# "Kij is the expected value of lnLKB in stratum i given sign j of the SAM during summer, 
# and f2 is the variance of lnLKB."
# So we test if there are difference in Ln(LKB) when SAM is + or -

LKB$logB <- log(LKB$biomass)
LKB$sam.sign <- ifelse(LKB$sam <0, "negative","positive")
LKB <- LKB %>% filter (LKB$gSSMU == "1" | LKB$gSSMU == "2")

two.way <- aov(logB ~ sam.sign + gSSMU, data = LKB)
summary(two.way)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# sam.sign      1   1.81   1.807   0.978 0.32433   
# gSSMU         1  19.44  19.440  10.526 0.00147 **
#   Residuals   141 260.40   1.847  
ggplot(LKB, aes(x = sam.sign, y = logB, col = gSSMU)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 2) +
  labs(x = "SAM annomaly", y = "Krill log(Biomass)") + 
  theme(axis.text = element_text(size = 12)) 

inter.way <- aov(logB ~ sam.sign*gSSMU, data = LKB)
summary(inter.way)
# Df Sum Sq Mean Sq F value Pr(>F)   
# sam.sign         1   1.81   1.807   0.975 0.3251   
# gSSMU            1  19.44  19.440  10.495 0.0015 **
#   sam.sign:gSSMU   1   1.07   1.074   0.580 0.4476   # Not significant
# Residuals      140 259.33   1.852 
boxplot(logB ~ sam.sign*gSSMU, data = LKB)



##########################################
# Relationship betwen LKB and Krill size ----------------------------------
##########################################

# Watters et al. 2020: "...we categorized estimates of LKB using a threshold of 1 ? 106 t (1 Mt)"
# But justification for this is just "Although it seems counterintuitive that the best case
# includes low LKB, some indices of penguin performance decrease when penguins forage on small krill..."
# We checked this using data frm the original paper by Hinke et al. 2007

krill.diet <- read.csv("./Supplementary Files/Krill_size_Hinke2007.csv") 

# First, we check if low LKB (1Mt) are correlated with smaller krill length in diet

# gSSMU1
Survey <- survey %>% filter(gSSMU == "1")
Krill <- krill.diet %>% filter(Location == "AB")
Survey.krill <- merge(x = Survey, y = Krill, by = "Year", all.x = TRUE)
Survey.krill$Year <- as.integer(Survey.krill$Year)
Survey.krill$LKB <- ifelse(Survey.krill$biomass > 1000000, ">1 Mt","< 1Mt")
Survey.krill <- Survey.krill %>% select(-gSSMU, -season, -Location)
LKB1.kr <- melt(Survey.krill, id = c("Year", "biomass", "LKB")) 
LKB1.kr<-na.omit(LKB1.kr)

ggplot(LKB1.kr, aes(x = LKB, y = value, col = variable)) +
  geom_boxplot(outlier.shape = NA) + 
  # geom_jitter(aes(col = variable, shape = variable), width = 0.5, size = 4) +
  geom_point(position = position_jitterdodge(), size = 4) +
  labs(x = "Local Krill Biomass", y = "Mean krill length in diet (mm)", title = "gSSMU1 summer") + 
  theme(axis.text = element_text(size = 12)) 

one.way <- aov(value ~ LKB, data = LKB1.kr)
summary(one.way)

two.way <- aov(value ~ LKB + variable, data = LKB1.kr)
summary(two.way)

inter.way <- aov(value ~ LKB*variable, data = LKB1.kr)
summary(inter.way)


# gSSMU2
Survey2 <- survey %>% filter(gSSMU == "2")
Krill2 <- krill.diet %>% filter(Location == "CS")
Survey2.krill <- merge(x = Survey2, y = Krill2, by = "Year", all.x = TRUE)
Survey2.krill$Year <- as.integer(Survey2.krill$Year)
Survey2.krill$LKB <- ifelse(Survey2.krill$biomass > 1000000, ">1 Mt","< 1Mt")
Survey2.krill <- Survey2.krill %>% select(-gSSMU, -season, -Location)
LKB2.kr <- melt(Survey2.krill, id = c("Year", "biomass", "LKB")) 
LKB2.kr<-na.omit(LKB2.kr)

ggplot(LKB2.kr, aes(x = LKB, y = value, col = variable)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 4) +
  labs(x = "Local Krill Biomass", y = "Mean krill size in diet (mm)", title = "gSSMU2 summer") + 
  theme(axis.text = element_text(size = 12)) 

one.way <- aov(value ~ LKB, data = LKB2.kr)
summary(one.way)

# For gSSMU2, small Krill is not related with low biomass (during summer) 

#All together
Survey <- survey %>% filter(gSSMU == "1" | gSSMU == "2")
Survey.krill <- merge(x = Survey, y = krill.diet, by = "Year", all.x = TRUE)
Survey.krill$Year <- as.integer(Survey.krill$Year)
Survey.krill$LKB <- ifelse(Survey.krill$biomass > 1000000, ">1 Mt","< 1Mt")
Survey.krill <- Survey.krill %>% select(-gSSMU, -season, -Location)
LKB.kr <- melt(Survey.krill, id = c("Year", "biomass", "LKB")) 
LKB.kr<-na.omit(LKB1.kr)

ggplot(LKB.kr, aes(x = LKB, y = value, col = variable)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), size = 4) +
  labs(x = "Local Krill Biomass", y = "Mean krill size in diet (mm)", title = "All summer data") + 
  theme(axis.text = element_text(size = 12))

# There is a trend for larger krill when LKB is low

inter.way <- aov(value ~ LKB*variable, data = LKB.kr)
summary(inter.way)


# Direct relationship between krill size and abundance
ggplot(LKB.kr, aes(x = biomass, y = value, col = variable, shape = variable)) +
  geom_point(size = 4) +
  labs(x = "Local Krill Biomass (ton)", y = "Mean krill size in diet (mm)", title = "All summer data") + 
  theme(axis.text = element_text(size = 12))
# No apparent correlation either

Survey$Year <- as.integer(Survey$Year)
ggplot(Survey, aes(x = Year, y = biomass, col = season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 4) +
  scale_y_continuous(breaks = seq(0, 10000000, 1000000)) +
  labs(y = "Local Krill Biomass (ton)", title = "All data gSSMU 1 & 2") + 
  theme(axis.text = element_text(size = 12))
# One problem is the bias influence of Winter surveys




# Penguin Indexes ---------------------------------------------------------

# Script from Watters et al
##############################
# generate the summer indices
##############################
# fledge weight (fwt)
# bigger indicates better summer
fwt <- S4
fwt$Year<-as.numeric(substr(fwt$YEAR,1,4))+1
fwt<-tapply(fwt$WT,list(fwt$Year,fwt$PROJECT,fwt$SPECIES),mean)
fwt<-data.frame(Year = rep(dimnames(fwt)[[1]],dim(fwt)[2]*dim(fwt)[3]),
                PROJECT=rep(rep(dimnames(fwt)[[2]],each=dim(fwt)[1]),dim(fwt)[3]),
                SPECIES=rep(dimnames(fwt)[[3]],each=dim(fwt)[1]*dim(fwt)[2]),
                fwt=c(fwt),stringsAsFactors = FALSE)
fwt$matchme<-paste(fwt$PROJECT,fwt$SPECIES,sep="|")
tt<-tapply(fwt$fwt,list(fwt$matchme),mean,na.rm=TRUE)
ttt<-tapply(fwt$fwt,list(fwt$matchme),sd,na.rm=TRUE)
mean.fwt<-tt[match(fwt$matchme,names(tt))]
sd.fwt<-ttt[match(fwt$matchme,names(ttt))]
fwt$std.mean.fwt<-(fwt$fwt-mean.fwt)/sd.fwt  # Why this index??
fwt<-fwt[,-c(5)]
# omits<-(fwt$SPECIES=="ADPE"&fwt$PROJECT=="CS")|(fwt$SPECIES=="CHPE"&fwt$PROJECT=="COPA")
# fwt<-fwt[!omits,]
fwt <- na.omit(fwt)
names(fwt)[5]<-"index1"
fwt$param = rep("FWT",dim(fwt)[1])
fwt$season = rep("S",dim(fwt)[1])
fwt$Year <- as.integer(fwt$Year)
# make stuff reference the correct "calendar year" for matching up with krill survey and catch data
# summer indices are relevant to the second year in the split-season designation
#print(str(fwt))


# post-hatch success (phs) (numbers of chicks creched/numbers of chicks hatched)
# bigger indicates better summer
phs <- S10
phs$phs<-phs$N_CRECHE/phs$N_CHICKS
phs$phs2<-log(phs$phs/(1-phs$phs))
phs$matchme<-paste(phs$PROJECT,phs$SPECIES,sep="|")
tt<-tapply(phs$phs2,list(phs$matchme),mean,na.rm=TRUE)
ttt<-tapply(phs$phs2,list(phs$matchme),sd,na.rm=TRUE)
mean.phs2 <- tt[match(phs$matchme,names(tt))]
sd.phs2 <- ttt[match(phs$matchme,names(ttt))]
phs$std.logit.phs2<-(phs$phs2-mean.phs2)/sd.phs2
phs <- phs[,-c(10)]
names(phs)[10]<-"index1"
phs$param=rep("PHS",dim(phs)[1])
phs$season=rep("S",dim(phs)[1])
# summer indices are relevant to the second year in the split-season designation
phs$Year<-as.numeric(substr(phs$YEAR,1,4))+1
#print(str(phs))


# trip duration (td)
# smaller indicates better summer (thus need to switch direction of index)
td <- S11
td<-td[,c(1:3,8)]
# next line is to make trip duration point in same direction as fwt and phs (max td is 59.95 for all trips)
# call this "revtd" for "reversed" trip duration
td[,4]<-60-td[,4]
names(td)[4]<-"revtd"
# summer indices are relevant to the second year in the split-season designation
td$Year<-as.numeric(substr(td$YEAR,1,4))+1
td<-tapply(td$revtd,list(td$Year,td$PROJECT,td$SPECIES),mean)
td<-data.frame(Year=rep(dimnames(td)[[1]],dim(td)[2]*dim(td)[3]),
               PROJECT=rep(rep(dimnames(td)[[2]],each=dim(td)[1]),dim(td)[3]),
               SPECIES=rep(dimnames(td)[[3]],each=dim(td)[1]*dim(td)[2]),
               revtd=c(td),stringsAsFactors = FALSE)
td$matchme<-paste(td$PROJECT,td$SPECIES,sep="|")
tt<-tapply(td$revtd,list(td$matchme),mean,na.rm=TRUE)
ttt<-tapply(td$revtd,list(td$matchme),sd,na.rm=TRUE)
mean.revtd<-tt[match(td$matchme,names(tt))]
sd.revtd<-ttt[match(td$matchme,names(ttt))]
td$std.revtd<-(td$revtd-mean.revtd)/sd.revtd
td<-td[,-c(5)]
names(td)[5]<-"index1"
#omits<-(td$SPECIES=="ADPE"&td$PROJECT=="CS")|(td$SPECIES=="CHPE"&td$PROJECT=="COPA")
#td<-td[!omits,]
td <- na.omit(td)
td$Year <- as.integer(td$Year)
td$param=rep("REVTD",dim(td)[1])
td$season=rep("S",dim(td)[1])
#print(str(td))

##############################
# generate the winter indices
##############################
#
# adult male mass at E1 lay (mml)
# bigger indicates better winter
ade1 <- S6
mml<-ade1[,c(1:3,5)]
# most winter indices (except rec) are relevant to the first year in the split-season designation
mml$Year<-as.numeric(substr(mml$YEAR,1,4))

mml<-tapply(mml$WT_MALE,list(mml$Year,mml$PROJECT,mml$SPECIES),mean,na.rm=TRUE)
mml<-data.frame(Year=rep(dimnames(mml)[[1]],dim(mml)[2]*dim(mml)[3]),
                PROJECT=rep(rep(dimnames(mml)[[2]],each=dim(mml)[1]),dim(mml)[3]),
                SPECIES=rep(dimnames(mml)[[3]],each=dim(mml)[1]*dim(mml)[2]),
                mml=c(mml),stringsAsFactors = FALSE)
mml$matchme<-paste(mml$PROJECT,mml$SPECIES,sep="|")
tt<-tapply(mml$mml,list(mml$matchme),mean,na.rm=TRUE)
ttt<-tapply(mml$mml,list(mml$matchme),sd,na.rm=TRUE)
mean.mml<-tt[match(mml$matchme,names(tt))]
sd.mml<-ttt[match(mml$matchme,names(ttt))]
mml$std.mean.mml<-(mml$mml-mean.mml)/sd.mml
mml<-mml[,-c(5)]
names(mml)[5]<-"index1"
#omits<-(mml$SPECIES=="ADPE"&mml$PROJECT=="CS")|(mml$SPECIES=="CHPE"&mml$PROJECT=="COPA")
#mml<-mml[!omits,]
mml <- na.omit(mml)
mml$Year <- as.integer(mml$Year)
mml$param=rep("MML",dim(mml)[1])
mml$season=rep("W",dim(mml)[1])
#print(str(mml))
#

# adult female mass at E1 lay (fml)
# bigger indicates better winter
ade1 <- S6
fml<-ade1[,c(1:3,6)]
# most winter indices (except rec) are relevant to the first year in the split-season designation
fml$Year<-as.numeric(substr(fml$YEAR,1,4))

fml<-tapply(fml$WT_FEMALE,list(fml$Year,fml$PROJECT,fml$SPECIES),mean,na.rm=TRUE)
fml<-data.frame(Year=rep(dimnames(fml)[[1]],dim(fml)[2]*dim(fml)[3]),
                PROJECT=rep(rep(dimnames(fml)[[2]],each=dim(fml)[1]),dim(fml)[3]),
                SPECIES=rep(dimnames(fml)[[3]],each=dim(fml)[1]*dim(fml)[2]),
                fml=c(fml),stringsAsFactors = FALSE)
fml$matchme<-paste(fml$PROJECT,fml$SPECIES,sep="|")
tt<-tapply(fml$fml,list(fml$matchme),mean,na.rm=TRUE)
ttt<-tapply(fml$fml,list(fml$matchme),sd,na.rm=TRUE)
mean.fml<-tt[match(fml$matchme,names(tt))]
sd.fml<-ttt[match(fml$matchme,names(ttt))]
fml$std.mean.fml<-(fml$fml-mean.fml)/sd.fml
fml<-fml[,-c(5)]
names(fml)[5]<-"index1"
#omits<-(fml$SPECIES=="ADPE"&fml$PROJECT=="CS")|(fml$SPECIES=="CHPE"&fml$PROJECT=="COPA")
#fml<-fml[!omits,]
fml <- na.omit(fml)
fml$Year <- as.integer(fml$Year)
fml$param=rep("fml",dim(fml)[1])
fml$season=rep("W",dim(fml)[1])
#print(str(mml))

# avg egg density using both eggs (egg)
# bigger indicates better winter
#function for scaling and centring the data
scale_this = function(x) {
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}
#create standardised centred index for egg density
  egg = S3 %>%
    dplyr::mutate(egg_den = c((E1_WT + E2_WT) / (E1_VOL + E2_VOL))) %>%
    dplyr::mutate(Year = as.numeric(str_sub(YEAR,1,4))) %>%
    drop_na() %>%
    dplyr::select(Year, PROJECT, SPECIES, egg_den) %>%
    dplyr::group_by(Year, PROJECT, SPECIES) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~
                               .x %>%
                               mutate(mean_egg = mean(egg_den)))) %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
    dplyr::mutate(match_me = paste(PROJECT,SPECIES,sep="|")) %>%
    dplyr::group_by(match_me) %>%
    nest() %>%
    dplyr::mutate(data = map(data, ~.x %>%
                               mutate(index1 = scale_this(mean_egg)))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(cols = c(data)) %>%
    dplyr::distinct(mean_egg, .keep_all=TRUE) %>%
    dplyr::mutate(param="EGG", season="W") %>%
    dplyr::select(Year, PROJECT, SPECIES, egg_den, index1, param, season)
####original EGG DENSITY code from Watters et al. (2020)----
  # avg egg density using both eggs (egg)
  # bigger indicates better winter
  e1e2 <- S3
  egg<-e1e2[,c(1:3)]
  egg$egg<-(e1e2[,5]+e1e2[,7])/(e1e2[,6]+e1e2[,8])
  # most winter indices (except rec) are relevant to the first year in the split-season designation
  egg$Year<-as.numeric(substr(egg$YEAR,1,4))
  egg<-tapply(egg$egg,list(egg$Year,egg$PROJECT,egg$SPECIES),mean,na.rm=TRUE)
  egg<-data.frame(Year=rep(dimnames(egg)[[1]],dim(egg)[2]*dim(egg)[3]),
                  PROJECT=rep(rep(dimnames(egg)[[2]],each=dim(egg)[1]),dim(egg)[3]),
                  SPECIES=rep(dimnames(egg)[[3]],each=dim(egg)[1]*dim(egg)[2]),
                  egg=c(egg),stringsAsFactors = FALSE)
  egg$matchme<-paste(egg$PROJECT,egg$SPECIES,sep="|")
  tt<-tapply(egg$egg,list(egg$matchme),mean,na.rm=TRUE)
  ttt<-tapply(egg$egg,list(egg$matchme),sd,na.rm=TRUE)
  mean.egg<-tt[match(egg$matchme,names(tt))]
  sd.egg<-ttt[match(egg$matchme,names(ttt))]
  egg$std.mean.egg<-(egg$egg-mean.egg)/sd.egg
  egg<-egg[,-c(5)]
  names(egg)[5]<-"index1"
  #omits<-(egg$SPECIES=="ADPE"&egg$PROJECT=="CS")|(egg$SPECIES=="CHPE"&egg$PROJECT=="COPA")
  egg <- na.omit(egg)
  egg$Year <- as.integer(egg$Year)
  egg$param=rep("EGG",dim(egg)[1])
  egg$season=rep("W",dim(egg)[1])
  #print(str(egg))
#### clutch initiation date (cid)----
# earlier indicates better winter
library(lubridate)

cid <- S2
# most winter indices (except rec) are relevant to the first year in the split-season designation
cid$Year<-as.numeric(substr(cid$YEAR,1,4))

# next line is to make CID point in same direction as other indices where bigger indicates better conditions (take diff from Dec 31)
# call this "revcid" for "reversed" CID
cid$revcid <- as.vector(as.POSIXlt(paste(substr(cid$YEAR,1,4),"-12-31",sep=""))-strptime(cid[,4],"%m/%e/%Y"))
cid$MEAN_CID <- as.Date(cid$MEAN_CID, tryFormats = c("%m/%d/%Y"))
cid$Julian <- lubridate::yday(cid$MEAN_CID)
cid$matchme<-paste(cid$PROJECT,cid$SPECIES,sep="|")
tt<-tapply(cid$revcid,list(cid$matchme),mean,na.rm=TRUE)
ttt<-tapply(cid$revcid,list(cid$matchme),sd,na.rm=TRUE)
mean.cid<-tt[match(cid$matchme,names(tt))]
sd.cid<-ttt[match(cid$matchme,names(ttt))]
cid$std.revcid<-(cid$revcid-mean.cid)/sd.cid
cid<-cid[,-c(1, 5,6,10)]
names(cid)[7]<-"index1"
cid$param=rep("REVCID",dim(cid)[1])
cid$season=rep("W",dim(cid)[1])
# uncomment next line if decide to remove gentoos because of their more plastic breeding phenology
#cid<-cid[cid$SPECIES!="GEPE",]
#print(str(cid))


# cohort recruitment (rec)
# bigger indicates better winter
rec <- S8[,c(1:5)]
#names(rec)[5]<-"rec"
rec$rec<-log(rec$RECRUITMENT/(1-rec$RECRUITMENT))
rec$matchme<-paste(rec$PROJECT,rec$SPECIES,sep="|")
tt<-tapply(rec$rec,list(rec$matchme),mean,na.rm=TRUE)
ttt<-tapply(rec$rec,list(rec$matchme),sd,na.rm=TRUE)
mean.rec<-tt[match(rec$matchme,names(tt))]
sd.rec<-ttt[match(rec$matchme,names(ttt))]
rec$std.logit.rec<-(rec$rec-mean.rec)/sd.rec
rec<-rec[,-c(7)]
names(rec)[7]<-"index1"
rec$param=rep("REC",dim(rec)[1])
rec$season=rep("W",dim(rec)[1])
# recruitment is relevant to the second year in the split-season designation (the winter of first independence)
#rec$cal.yr<-as.numeric(substr(rec$YEAR,1,4))+1
names(rec)[1] <- "Year"
#print(str(rec))



# Relationship KRILL size and Penguin Indexes -----------------------------

#Krill size is only summer data; varies depending on Species and Location
#krill.diet <- read.csv("Krill_size_Hinke2007.csv") 
krill.diet$Location <- str_replace_all(krill.diet$Location, c(AB = "COPA"))
names(krill.diet) <- c("Year", "ADPE", "CHPE", "GEPE", "PROJECT")
krill.spp <- melt(krill.diet, id = c("Year", "PROJECT"))
names(krill.spp) <- c("Year", "PROJECT", "SPECIES", "Krill.size")
krill.spp <- na.omit(krill.spp)

fwt.k <- merge(x = fwt, y = krill.spp, by = c("Year", "PROJECT", "SPECIES"))
phs.k <- merge(x = phs, y = krill.spp, by = c("Year", "PROJECT", "SPECIES"))
td.k <- merge(x = td, y = krill.spp, by = c("Year", "PROJECT", "SPECIES"))

fwt2 <- fwt.k %>% select(Year, PROJECT, SPECIES, index1, param, Krill.size)
phs2 <- phs.k %>% select(Year, PROJECT, SPECIES, index1, param, Krill.size)
td2 <- td.k %>% select(Year, PROJECT, SPECIES, index1, param, Krill.size)

Peng.krill <- rbind(fwt2, phs2, td2)

CS <- Peng.krill %>% filter(PROJECT == "CS")
AB <- Peng.krill %>% filter(PROJECT == "COPA")

#CS
ggplot(CS, aes(x = Krill.size, y = index1, col = SPECIES, shape = SPECIES)) +
  geom_point(size = 4) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(x = "Mean krill size (mm)", y = "Summer Indexes", title = "CAPE SHIRREFF") + 
  theme(axis.text = element_text(size = 12))

CS.lm <- lm(Krill.size ~ index1:SPECIES, data = CS)
summary(CS.lm)

two.way1 <- aov(index1 ~ oni.class + param, data = Peng.summer)
summary(two.way)

#AB
ggplot(AB, aes(x = Krill.size, y = index1, col = SPECIES, shape = SPECIES)) +
  geom_point(size = 4) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(x = "Mean krill size (mm)", y = "Summer Indexes", title = "ADMIRALTY BAY") + 
  theme(axis.text = element_text(size = 12))

ab.lm <- lm(Krill.size ~ index1:SPECIES, data = AB)
summary(ab.lm)

two.way1 <- aov(index1 ~ oni.class + param, data = Peng.summer)
summary(two.way)


# Relationship ONI and Penguin Indexes ------------------------------------

#Estimate Summer and Winter Indexes - ONI
oni <- S7
oni$yr<-ifelse(is.element(oni$SEAS,c("OND","NDJ")),oni$YR+1,oni$YR)
oni$season<-ifelse(is.element(oni$SEAS,c("OND","NDJ","DJF","JFM")),"S",NA)
oni$season<-ifelse(is.element(oni$SEAS,c("AMJ","MJJ","JJA","JAS")),"W",oni$season)
oni<-na.omit(oni)
oni<-tapply(oni$ANOM,list(oni$yr,oni$season),mean)
oni<-data.frame(yr=rep(dimnames(oni)[[1]],2),season=rep(dimnames(oni)[[2]],each=dim(oni)[1]),oni=c(oni))
names(oni)<-c("Year","season","oni")
oni$Year <- as.integer(oni$Year)
oni$oni.class <- ifelse(oni$oni <= -0.5, "Cool","Neutral")
oni$oni.class<-ifelse(oni$oni >= 0.5, "Warm", oni$oni.class)


# Summer Indexes
oni.s <- oni %>% filter(season == "S")
fwt.oni <- merge(x = fwt, y = oni.s, by = "Year", all.x = TRUE)
phs.oni <- merge(x = phs, y = oni.s, by = "Year", all.x = TRUE)
td.oni <- merge(x = td, y = oni.s, by = "Year", all.x = TRUE)
fwt2 <- fwt.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
phs2 <- phs.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
td2 <- td.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)

Peng.summer <- rbind(fwt2, phs2, td2)
ggplot(Peng.summer, aes(x = oni.class, y = index1, col = param)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 4) +
  labs(x = "ONI", y = "Summer Indexes", title = "Summer Indexes: Parameters") + 
  theme(axis.text = element_text(size = 12))

two.way1 <- aov(index1 ~ oni.class + param, data = Peng.summer)
summary(two.way)

tukey.two.way<-TukeyHSD(two.way1); tukey.two.way

inter.way1 <- aov(index1 ~ oni.class*param, data = Peng.summer)
summary(inter.way)


ggplot(Peng.summer, aes(x = oni.class, y = index1, col = SPECIES)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 4) +
  labs(x = "ONI", y = "Summer Indexes", title = "Summer Indexes: Species") + 
  theme(axis.text = element_text(size = 12))

two.way2 <- aov(index1 ~ oni.class + SPECIES, data = Peng.summer)
summary(two.way)

tukey.two.way<-TukeyHSD(two.way); tukey.two.way

inter.way2 <- aov(index1 ~ oni.class*SPECIES, data = Peng.summer)
summary(inter.way)
tukey.inter<-TukeyHSD(inter.way2); tukey.inter

model.set <- list(two.way1, two.way2, inter.way1, inter.way2)
model.names <- c("two.way param", "two.way Species", "Inter.param", "Inter.Species")
aictab(model.set, modnames = model.names)


# Winter Indexes
oni.w <- oni %>% filter(season == "W")
mml.oni <- merge(x = mml, y = oni.w, by = "Year", all.x = TRUE)
mml2 <- mml.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
fml.oni <- merge(x = fml, y = oni.w, by = "Year", all.x = TRUE)
fml2 <- fml.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
egg.oni <- merge(x = egg, y = oni.w, by = "Year", all.x = TRUE)
egg2 <- egg.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
cid.oni <- merge(x = cid, y = oni.w, by = "Year", all.x = TRUE)
cid2 <- cid.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)
rec.oni <- merge(x = rec, y = oni.w, by = "Year", all.x = TRUE)
rec2 <- rec.oni %>% select(Year, PROJECT, SPECIES, index1, param, oni.class)

Peng.winter <- rbind(mml2, fml2, egg2, cid2, rec2)

#by parameters
ggplot(Peng.winter, aes(x = oni.class, y = index1, col = param)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 2) +
  labs(x = "ONI", y = "Winter Indexes", title = "Winter Indexes: Parameters") + 
  theme(axis.text = element_text(size = 12))

two.way1 <- aov(index1 ~ oni.class + param, data = Peng.winter)
summary(two.way1)

tukey.two.way<-TukeyHSD(two.way1); tukey.two.way

inter.way1 <- aov(index1 ~ oni.class*param, data = Peng.winter)
summary(inter.way1)
tukey.inter<-TukeyHSD(inter.way1); tukey.inter


#by species
ggplot(Peng.winter, aes(x = oni.class, y = index1, col = SPECIES)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 2) +
  labs(x = "ONI", y = "Winter Indexes", title = "Winter Indexes: Species") + 
  theme(axis.text = element_text(size = 12))

two.way2 <- aov(index1 ~ oni.class + SPECIES, data = Peng.winter)
summary(two.way2)

tukey.two.way<-TukeyHSD(two.way2); tukey.two.way

inter.way2 <- aov(index1 ~ oni.class*SPECIES, data = Peng.winter)
summary(inter.way2)
tukey.inter<-TukeyHSD(inter.way2); tukey.inter

model.set <- list(two.way1, two.way2, inter.way1, inter.way2)
model.names <- c("two.way param", "two.way Species", "Inter.param", "Inter.Species")
aictab(model.set, modnames = model.names)
