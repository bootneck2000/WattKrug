#### MODEL EXECUTION
## WATTERS ET AL (2020)

## Conditions: 
# A. Survey data filtered by transect (nm) coverage 
# B. Not imputed summer LKB
# C. Penguin distribution adjusted to specific SSMUs

library(tidyverse)
library(lattice)
library(dplyr)
select <- dplyr::select
filter <- dplyr::filter

# Jump to Step 3 if already ran steps 1 & 2

#### First Step, load variables & Update 'survey' data ####
read_dir <- "./Supplementary Files/"

make.localhr.data<-function(trim=1,plot.winter=FALSE){
  ###########################################################################################################
  # generate the summer indices
  ###########################################################################################################
  #
  # fledge weight (fwt)
  # bigger indicates better summer
  fwt<-read.csv(file.path(read_dir,"fweight.csv"),header=TRUE,stringsAsFactors = FALSE)
  fwt<-tapply(fwt$WT,list(fwt$YEAR,fwt$PROJECT,fwt$SPECIES),mean)
  fwt<-data.frame(YEAR=rep(dimnames(fwt)[[1]],dim(fwt)[2]*dim(fwt)[3]),
                  PROJECT=rep(rep(dimnames(fwt)[[2]],each=dim(fwt)[1]),dim(fwt)[3]),
                  SPECIES=rep(dimnames(fwt)[[3]],each=dim(fwt)[1]*dim(fwt)[2]),
                  fwt=c(fwt),stringsAsFactors = FALSE)
  fwt$matchme<-paste(fwt$PROJECT,fwt$SPECIES,sep="|")
  tt<-tapply(fwt$fwt,list(fwt$matchme),mean,na.rm=TRUE)
  ttt<-tapply(fwt$fwt,list(fwt$matchme),sd,na.rm=TRUE)
  mean.fwt<-tt[match(fwt$matchme,names(tt))]
  sd.fwt<-ttt[match(fwt$matchme,names(ttt))]
  fwt$std.mean.fwt<-(fwt$fwt-mean.fwt)/sd.fwt
  fwt<-fwt[,-c(4:5)]
  #omits<-(fwt$SPECIES=="ADPE"&fwt$PROJECT=="CS")|(fwt$SPECIES=="CHPE"&fwt$PROJECT=="COPA")
  #fwt<-fwt[!omits,]
  names(fwt)[4]<-"index"
  fwt$param=rep("FWT",dim(fwt)[1])
  fwt$season=rep("S",dim(fwt)[1])
  # make stuff reference the correct "calendar year" for matching up with krill survey and catch data
  # summer indices are relevant to the second year in the split-season designation
  fwt$cal.yr<-as.numeric(substr(fwt$YEAR,1,4))+1
  #print(str(fwt))
  #
  ###########################################################################################################
  # post-hatch success (phs) (numbers of chicks creched/numbers of chicks hatched)
  # bigger indicates better summer
  phs<-read.csv(file.path(read_dir,"success.csv"),header=TRUE,stringsAsFactors = FALSE)
  phs$phs<-phs$N_CRECHE/phs$N_CHICKS
  phs$phs<-log(phs$phs/(1-phs$phs))
  phs$matchme<-paste(phs$PROJECT,phs$SPECIES,sep="|")
  tt<-tapply(phs$phs,list(phs$matchme),mean,na.rm=TRUE)
  ttt<-tapply(phs$phs,list(phs$matchme),sd,na.rm=TRUE)
  mean.phs<-tt[match(phs$matchme,names(tt))]
  sd.phs<-ttt[match(phs$matchme,names(ttt))]
  phs$std.logit.phs<-(phs$phs-mean.phs)/sd.phs
  phs<-phs[,-c(4:9)]
  names(phs)[4]<-"index"
  phs$param=rep("PHS",dim(phs)[1])
  phs$season=rep("S",dim(phs)[1])
  # summer indices are relevant to the second year in the split-season designation
  phs$cal.yr<-as.numeric(substr(phs$YEAR,1,4))+1
  #print(str(phs))
  #
  ###########################################################################################################
  # trip duration (td)
  # smaller indicates better summer (thus need to switch direction of index)
  td<-read.csv(file.path(read_dir,"tripduration.csv"),header=TRUE,stringsAsFactors = FALSE)
  td<-td[,c(1:3,8)]
  # next line is to make trip duration point in same direction as fwt and phs (max td is 59.95 for all trips)
  # call this "revtd" for "reversed" trip duration
  td[,4]<-60-td[,4]
  names(td)[4]<-"revtd"
  td<-tapply(td$revtd,list(td$YEAR,td$PROJECT,td$SPECIES),mean)
  td<-data.frame(YEAR=rep(dimnames(td)[[1]],dim(td)[2]*dim(td)[3]),
                 PROJECT=rep(rep(dimnames(td)[[2]],each=dim(td)[1]),dim(td)[3]),
                 SPECIES=rep(dimnames(td)[[3]],each=dim(td)[1]*dim(td)[2]),
                 revtd=c(td),stringsAsFactors = FALSE)
  td$matchme<-paste(td$PROJECT,td$SPECIES,sep="|")
  tt<-tapply(td$revtd,list(td$matchme),mean,na.rm=TRUE)
  ttt<-tapply(td$revtd,list(td$matchme),sd,na.rm=TRUE)
  mean.revtd<-tt[match(td$matchme,names(tt))]
  sd.revtd<-ttt[match(td$matchme,names(ttt))]
  td$std.revtd<-(td$revtd-mean.revtd)/sd.revtd
  td<-td[,-c(4:5)]
  names(td)[4]<-"index"
  #omits<-(td$SPECIES=="ADPE"&td$PROJECT=="CS")|(td$SPECIES=="CHPE"&td$PROJECT=="COPA")
  #td<-td[!omits,]
  td$param=rep("REVTD",dim(td)[1])
  td$season=rep("S",dim(td)[1])
  # summer indices are relevant to the second year in the split-season designation
  td$cal.yr<-as.numeric(substr(td$YEAR,1,4))+1
  #print(str(td))
  #
  ###########################################################################################################
  # generate the winter indices
  ###########################################################################################################
  #
  # adult male mass at E1 lay (mml)
  # bigger indicates better winter
  ade1<-read.csv(file.path(read_dir,"massatlay.csv"),header=TRUE,stringsAsFactors = FALSE)
  mml<-ade1[,c(1:3,5)]
  mml<-tapply(mml$WT_MALE,list(mml$YEAR,mml$PROJECT,mml$SPECIES),mean,na.rm=TRUE)
  mml<-data.frame(YEAR=rep(dimnames(mml)[[1]],dim(mml)[2]*dim(mml)[3]),
                  PROJECT=rep(rep(dimnames(mml)[[2]],each=dim(mml)[1]),dim(mml)[3]),
                  SPECIES=rep(dimnames(mml)[[3]],each=dim(mml)[1]*dim(mml)[2]),
                  mml=c(mml),stringsAsFactors = FALSE)
  mml$matchme<-paste(mml$PROJECT,mml$SPECIES,sep="|")
  tt<-tapply(mml$mml,list(mml$matchme),mean,na.rm=TRUE)
  ttt<-tapply(mml$mml,list(mml$matchme),sd,na.rm=TRUE)
  mean.mml<-tt[match(mml$matchme,names(tt))]
  sd.mml<-ttt[match(mml$matchme,names(ttt))]
  mml$std.mean.mml<-(mml$mml-mean.mml)/sd.mml
  mml<-mml[,-c(4:5)]
  names(mml)[4]<-"index"
  #omits<-(mml$SPECIES=="ADPE"&mml$PROJECT=="CS")|(mml$SPECIES=="CHPE"&mml$PROJECT=="COPA")
  #mml<-mml[!omits,]
  mml$param=rep("MML",dim(mml)[1])
  mml$season=rep("W",dim(mml)[1])
  # most winter indices (except rec) are relevant to the first year in the split-season designation
  mml$cal.yr<-as.numeric(substr(mml$YEAR,1,4))
  #print(str(mml))
  #
  ###########################################################################################################
  #
  # adult female mass at E1 lay (fml)
  # bigger indicates better winter
  fml<-ade1[,c(1:3,6)]
  fml<-tapply(fml$WT_FEMALE,list(fml$YEAR,fml$PROJECT,fml$SPECIES),mean,na.rm=TRUE)
  fml<-data.frame(YEAR=rep(dimnames(fml)[[1]],dim(fml)[2]*dim(fml)[3]),
                  PROJECT=rep(rep(dimnames(fml)[[2]],each=dim(fml)[1]),dim(fml)[3]),
                  SPECIES=rep(dimnames(fml)[[3]],each=dim(fml)[1]*dim(fml)[2]),
                  fml=c(fml),stringsAsFactors = FALSE)
  fml$matchme<-paste(fml$PROJECT,fml$SPECIES,sep="|")
  tt<-tapply(fml$fml,list(fml$matchme),mean,na.rm=TRUE)
  ttt<-tapply(fml$fml,list(fml$matchme),sd,na.rm=TRUE)
  mean.fml<-tt[match(fml$matchme,names(tt))]
  sd.fml<-ttt[match(fml$matchme,names(ttt))]
  fml$std.mean.fml<-(fml$fml-mean.fml)/sd.fml
  fml<-fml[,-c(4:5)]
  names(fml)[4]<-"index"
  #omits<-(fml$SPECIES=="ADPE"&fml$PROJECT=="CS")|(fml$SPECIES=="CHPE"&fml$PROJECT=="COPA")
  #fml<-fml[!omits,]
  fml$param=rep("FML",dim(fml)[1])
  fml$season=rep("W",dim(fml)[1])
  # most winter indices (except rec) are relevant to the first year in the split-season designation
  fml$cal.yr<-as.numeric(substr(fml$YEAR,1,4))
  #print(str(fml))
  #
  ###########################################################################################################
  #
  # avg egg density using both eggs (egg)
  # bigger indicates better winter
  e1e2<-read.csv(file.path(read_dir,"egg.csv"),header=TRUE,stringsAsFactors = FALSE)
  egg<-e1e2[,c(1:3)]
  egg$egg<-(e1e2[,5]+e1e2[,7])/(e1e2[,6]+e1e2[,8])
  egg<-tapply(egg$egg,list(egg$YEAR,egg$PROJECT,egg$SPECIES),mean,na.rm=TRUE)
  egg<-data.frame(YEAR=rep(dimnames(egg)[[1]],dim(egg)[2]*dim(egg)[3]),
                  PROJECT=rep(rep(dimnames(egg)[[2]],each=dim(egg)[1]),dim(egg)[3]),
                  SPECIES=rep(dimnames(egg)[[3]],each=dim(egg)[1]*dim(egg)[2]),
                  egg=c(egg),stringsAsFactors = FALSE)
  egg$matchme<-paste(egg$PROJECT,egg$SPECIES,sep="|")
  tt<-tapply(egg$egg,list(egg$matchme),mean,na.rm=TRUE)
  ttt<-tapply(egg$egg,list(egg$matchme),sd,na.rm=TRUE)
  mean.egg<-tt[match(egg$matchme,names(tt))]
  sd.egg<-ttt[match(egg$matchme,names(ttt))]
  egg$std.mean.egg<-(egg$egg-mean.egg)/sd.egg
  egg<-egg[,-c(4:5)]
  names(egg)[4]<-"index"
  #omits<-(egg$SPECIES=="ADPE"&egg$PROJECT=="CS")|(egg$SPECIES=="CHPE"&egg$PROJECT=="COPA")
  #egg<-egg[!omits,]
  egg$param=rep("EGG",dim(egg)[1])
  egg$season=rep("W",dim(egg)[1])
  # most winter indices (except rec) are relevant to the first year in the split-season designation
  egg$cal.yr<-as.numeric(substr(egg$YEAR,1,4))
  #print(str(egg))
  #
  ###########################################################################################################
  #
  # clutch initiation date (cid)
  # earlier indicates better winter
  cid<-read.csv(file.path(read_dir,"cid.csv"),header=TRUE,stringsAsFactors = FALSE)[,1:4]
  # next line is to make CID point in same direction as other indices where bigger indicates better conditions (take diff from Dec 31)
  # call this "revcid" for "reversed" CID
  cid[,4]<-as.vector(as.POSIXlt(paste(substr(cid$YEAR,1,4),"-12-31",sep=""))-strptime(cid[,4],"%m/%e/%Y"))
  names(cid)[4]<-"revcid"
  cid$matchme<-paste(cid$PROJECT,cid$SPECIES,sep="|")
  tt<-tapply(cid$revcid,list(cid$matchme),mean,na.rm=TRUE)
  ttt<-tapply(cid$revcid,list(cid$matchme),sd,na.rm=TRUE)
  mean.cid<-tt[match(cid$matchme,names(tt))]
  sd.cid<-ttt[match(cid$matchme,names(ttt))]
  cid$std.revcid<-(cid$revcid-mean.cid)/sd.cid
  cid<-cid[,-c(4:5)]
  names(cid)[4]<-"index"
  cid$param=rep("REVCID",dim(cid)[1])
  cid$season=rep("W",dim(cid)[1])
  # most winter indices (except rec) are relevant to the first year in the split-season designation
  cid$cal.yr<-as.numeric(substr(cid$YEAR,1,4))
  # uncomment next line if decide to remove gentoos because of their more plastic breeding phenology
  #cid<-cid[cid$SPECIES!="GEPE",]
  #print(str(cid))
  #
  ###########################################################################################################
  #
  # cohort recruitment (rec)
  # bigger indicates better winter
  rec<-read.csv(file.path(read_dir,"recruitment.csv"),header=TRUE,stringsAsFactors = FALSE)[,1:4]
  names(rec)[4]<-"rec"
  rec$rec<-log(rec$rec/(1-rec$rec))
  rec$matchme<-paste(rec$PROJECT,rec$SPECIES,sep="|")
  tt<-tapply(rec$rec,list(rec$matchme),mean,na.rm=TRUE)
  ttt<-tapply(rec$rec,list(rec$matchme),sd,na.rm=TRUE)
  mean.rec<-tt[match(rec$matchme,names(tt))]
  sd.rec<-ttt[match(rec$matchme,names(ttt))]
  rec$std.logit.rec<-(rec$rec-mean.rec)/sd.rec
  rec<-rec[,-c(4:5)]
  names(rec)[4]<-"index"
  rec$param=rep("REC",dim(rec)[1])
  rec$season=rep("W",dim(rec)[1])
  # recruitment is relevant to the second year in the split-season designation (the winter of first independence)
  rec$cal.yr<-as.numeric(substr(rec$YEAR,1,4))+1
  #print(str(rec))
  #
  ###########################################################################################################
  # read in the krill survey and fishery data
  ###########################################################################################################
  #
  # krill survey biomass
  survey<-read.csv(file.path(read_dir,"krillsurveywithJoinville.csv"),header=TRUE,stringsAsFactors = FALSE)
  # use next line if want to filter acoustic data to have minimum number of miles (comment out if not desired)
  # as per CSR, 80 nmi would be about equivalent of 2 tracklines in the Bransfield
  #survey<-survey[survey$nmi.count>=80,]
  # could try changing "biomass" in following line to "mean.density.gm2" or "median.density.gm2" but haven't done that

  ### Introduced Change: 
  # Filter low transect coverage (<10th percentil)
  survey_filtered <- survey %>%
    filter(gSSMU %in% c(1, 2))
  percentiles <- survey_filtered %>%
    group_by(gSSMU) %>%
    summarise(p10 = quantile(nmi.count, probs = 0.10, na.rm = TRUE))
  percentiles$p10 <- round(percentiles$p10,0)
  survey_low_nmi <- survey_filtered %>%
    left_join(percentiles, by = "gSSMU") %>%
    filter(nmi.count >= p10)
  survey <- survey_low_nmi
  rm(percentiles, survey_filtered, survey_low_nmi)
  ### END changes
  
  ### New change: Estimate SSMU-LKB
  
  # Read SSMU table
  SSMUs<-read.csv(file.path(read_dir,"SSMU48_area.csv"),header=TRUE,stringsAsFactors = FALSE)
  
  survey.mod <- survey %>%
    left_join(SSMUs %>% select(gSSMU, SSMU, Area.fraction), by = "gSSMU",
              relationship = "many-to-many") %>%
    mutate(biomass2 = biomass * Area.fraction)
  
  survey <- survey.mod %>%
    select(-biomass) %>%
    rename(biomass = biomass2)
    ### END changes
  
  survey<-tapply(survey$biomass,list(survey$Year,survey$SSMU),mean,na.rm=TRUE) # change gSSMU by SSMU
  survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                     SSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                     survey=c(survey),stringsAsFactors = FALSE)
  survey$season<-ifelse(survey$cal.yr<2012,"S","W")
  # use next line if want to remove winter survey data altogether (comment out if not desired)
  #survey<-survey[survey$season=="S",]
  survey$matchme<-paste(survey$cal.yr,survey$season,survey$SSMU,sep="|")
  #print(str(survey))
  rm(survey.mod)
  #
  ###########################################################################################################
  #
  # krill fishery catches
  fishery<-read.csv(file.path(read_dir,"c1.csv"),header=TRUE,stringsAsFactors = FALSE)
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

  fishery<- fishery %>% rename(SSMU = AssignedSSMU) # change column name to SSMU
  
  fishery<-tapply(fishery$TotalCatch,list(fishery$CalendarYear,fishery$SSMU,fishery$season),sum)
  fishery<-data.frame(cal.yr=rep(dimnames(fishery)[[1]],dim(fishery)[2]*dim(fishery)[3]),
                      SSMU=rep(rep(dimnames(fishery)[[2]],each=dim(fishery)[1]),dim(fishery)[3]),
                      season=rep(dimnames(fishery)[[3]],each=dim(fishery)[1]*dim(fishery)[2]),
                      catch=c(fishery),stringsAsFactors = FALSE)
  fishery$cal.yr<-as.numeric(as.character(fishery$cal.yr))
  fishery$matchme<-paste(fishery$cal.yr,fishery$season,fishery$SSMU,sep="|")
  
  fishery$catch[is.na(fishery$catch)] <- 0 # Added Change
  fishery$cal.yr <- as.character(fishery$cal.yr)
  
  #print(str(fishery))
  #
  ###########################################################################################################
  # now match predator data with krill data
  out<-rbind(fwt,phs,td,mml,fml,egg,cid,rec,make.row.names=FALSE)
  # all birds from Copa always forage in gSSMU 1 (Bransfield SSMUs)
  # CHPE from Cape Shirreff always forage in gSSMU 2 (Drake Passage SSMUs)
  # GEPE from Cape Shirreff forage in gSSMU 2 during summer and gSSMU 1 during winter
  #out$gSSMU<-ifelse(out$PROJECT=="COPA",1,
  #                  ifelse(out$SPECIES=="CHPE",2,
  #                         ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="S",2,1)))
  
  ### Changes in Spatial Distribution
  # Summer: COPA penguins only feed in APBSE
  # Summer: SHIRREFF penguins only feed in APDPW
  # Winter: COPA: Chinstraps feed in APDPW, APDPE, APEI
  # Winter: COPA: Adelies feed in APBSE, APDPE, APEI
  # Winter: COPA: Gentoos feed in APBSE
  # Winter: SHIRREFF: Chinstraps feed in APDPW, APDPE, APEI
  # Winter: SHIRREFF: Gentoos feed in APBSW
  
  out<-out[!is.na(out$index),] # remove empty index data
  
  out$SSMU<-rep(NA,dim(out)[1])
  out$SSMU<-ifelse(out$SPECIES=="ADPE"&out$PROJECT=="COPA"&out$season=="S","APBSE",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="ADPE"&out$PROJECT=="COPA"&out$season=="W","APBSE|APDPE|APEI",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="S","APBSE",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="W","APDPW|APDPE|APEI",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="S","APBSE",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="W","APBSE",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="S","APDPW",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="W","APDPW|APDPE|APEI",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="S","APDPW",out$SSMU)
  out$SSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="W","APBSW",out$SSMU)
  #
  out$matchme<-paste(out$cal.yr,out$season,out$SSMU,sep="|")
  # out$survey<-survey$survey[match(out$matchme,survey$matchme)]
  # out$catch<-fishery$catch[match(out$matchme,fishery$matchme)]
  
  ### Change code to handle multiple SSMU
  # First, split the SSMU column in 'out' to handle multiple SSMUs
  out$cal.yr <- as.character(out$cal.yr)
  
  out_expanded <- out %>%
    mutate(SSMU_list = strsplit(as.character(SSMU), "\\|")) %>%
    unnest(SSMU_list) %>%
    # Join with survey data
    left_join(survey %>% select(SSMU, cal.yr, season, survey), 
              by = c("SSMU_list" = "SSMU", "cal.yr", "season")) %>%
    # Join with fishery data
    left_join(fishery %>% select(SSMU, cal.yr, season, catch), 
              by = c("SSMU_list" = "SSMU", "cal.yr", "season")) %>%
    # Sum survey and catch values for each original row
    group_by(across(-c(SSMU_list, survey, catch))) %>%
    summarise(
      survey = sum(survey, na.rm = TRUE),
      catch = sum(catch, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Replace 'out' with the result
  out <- out_expanded  

  #
  # out<-out[!is.na(out$SSMU),]
  #
  ###########################################################################################################
  # pull in the environmental indices
  ###########################################################################################################
  #
  # SOUTHERN ANNULAR MODE
  sam<-read.csv(file.path(read_dir,"sam.csv"))
  names(sam)<-c("yr","mo","sam")
  sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
  sam$YEAR<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
  sam<-tapply(sam$sam,list(sam$YEAR,sam$season),mean)
  sam<-data.frame(YEAR=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))
  out$sam<-sam$sam[match(paste(out$cal.yr,out$season,sep="|"),paste(sam$YEAR,sam$season,sep="|"))]
  out$sam.sign<-ifelse(out$sam<0,"Neg","Pos")
  #
  # OCEANIC NINO INDEX
  oni<-read.csv(file.path(read_dir,"oni.csv"),stringsAsFactors = FALSE)
  oni$yr<-ifelse(is.element(oni$SEAS,c("OND","NDJ")),oni$YR+1,oni$YR)
  oni$season<-ifelse(is.element(oni$SEAS,c("OND","NDJ","DJF","JFM")),"S",NA)
  oni$season<-ifelse(is.element(oni$SEAS,c("AMJ","MJJ","JJA","JAS")),"W",oni$season)
  oni<-na.omit(oni)
  oni<-tapply(oni$ANOM,list(oni$yr,oni$season),mean)
  oni<-data.frame(yr=rep(dimnames(oni)[[1]],2),season=rep(dimnames(oni)[[2]],each=dim(oni)[1]),oni=c(oni))
  out$oni<-oni$oni[match(paste(out$cal.yr,out$season,sep="|"),paste(oni$yr,oni$season,sep="|"))]
  out$oni.class<-ifelse(out$oni <= -0.5, "Cool","Neutral")
  out$oni.class<-ifelse(out$oni >=0.5, "Warm",out$oni.class)
  #
  #
  # some clean up
  #
  out$catch[is.na(out$catch)]<-0
  out<-out[!is.na(out$sam),]
  out<-out[!is.nan(out$index),]
  out<-out[!is.na(out$index),]
  # will not try to impute missing winter surveys
  # but will keep winter performance indices if want to plot them
  if(!plot.winter){out<-out[!(is.na(out$survey)&out$season=="W"),]}
  #
  # if require minimum number of data points per study
  if(!is.null(trim)){
    study<-as.numeric(factor(paste(out$PROJECT,out$SPECIES,out$param,sep="|")))
    study.n<-table(study)
    keepers<-as.numeric(as.vector(dimnames(study.n[study.n>trim])[[1]]))
    out<-out[is.element(study,keepers),]
  }
  out
}

junk<-make.localhr.data()

junk<-filter(junk,survey>0)
junk$bclass<-ifelse(junk$survey<=1000000,1,2)
junk$hrclass<-ifelse(junk$catch/junk$survey<=0.01,1,ifelse(junk$catch/junk$survey>=0.1,3,2))
junk$oniclass<-ifelse(junk$oni.class=="Cool",1,ifelse(junk$oni.class=="Warm",3,2))



#### Second - Run Model ####

modelstring<-"
  # George Watters -- April 2019

  # oniclass are environment classes based on ONI where
  # 1 = cool
  # 2 = neutral
  # 3 = warm

  model{

    #####################################################################
    #
    # model for penguin performance
    #
    #####################################################################


    # the design matrix with sum-to-zero constraints
    for(i in 1:nobs){
      X[i,1]<-1.0     # intercept
      X[i,2]<-equals(bclass[i],2)-equals(bclass[i],1) # b2
      X[i,3]<-equals(hrclass[i],2)-equals(hrclass[i],1) # hr2
      X[i,4]<-equals(hrclass[i],3)-equals(hrclass[i],1) # hr3
      X[i,5]<-equals(oniclass[i],2)-equals(oniclass[i],1) # o2
      X[i,6]<-equals(oniclass[i],3)-equals(oniclass[i],1) # o3
    }
  

    # the likelihood
    for(i in 1:nobs){
      index[i]~dnorm(mu.index[i],tau.index)
      mu.index[i] <- inprod(X[i,],beta[])
    }
    

    # priors

    beta[1]~dnorm(0, 0.0001)
    
    beta[2]~dnorm(0, 0.0001)
    
    beta[3]~dnorm(0, 0.0001)
    beta[4]~dnorm(0, 0.0001)

    beta[5]~dnorm(0, 0.0001)
    beta[6]~dnorm(0, 0.0001)
    

    # half-cauchy for variation among indices
    tau.index<-pow(sd.index,-2)
    #sd.index~dunif(0,10)
    sd.index~dt(0,t.tau.index,1)T(0,)
    t.tau.index<-pow(t.sd.index,-2)
    # hyperprior for half-cauchy scale
    t.sd.index~dunif(0,2)


    # derived quantities
    # first the design matrix for easily interpreting effects
    # row 1 -- ONI=cool, LKB<=1Mt, LHR<=0.01 (reference case)
    # row 2 -- ONI=cool, LKB>1Mt, 0.01<LHR<0.1
    # row 3 -- ONI=cool, LKB<=1Mt, LHR>=0.1
    # row 4 -- ONI=cool, LKB>1Mt, LHR<=0.01
    # row 5 -- ONI=cool, LKB<=1Mt, 0.01<LHR<0.1
    # row 6 -- ONI=cool, LKB>1Mt, LHR>=0.1
    # row 7 -- ONI=neutral, LKB<=1Mt, LHR<=0.01
    # row 8 -- ONI=neutral, LKB>1Mt, 0.01<LHR<0.1
    # row 9 -- ONI=neutral, LKB<=1Mt, LHR>=0.1
    # row 10 -- ONI=neutral, LKB>1Mt, LHR<=0.01
    # row 11 -- ONI=neutral, LKB<=1Mt, 0.01<LHR<0.1
    # row 12 -- ONI=neutral, LKB>1Mt, LHR>=0.1 (worst case)
    # row 13 -- ONI=warm, LKB<=1Mt, LHR<=0.01
    # row 14 -- ONI=warm, LKB>1Mt, 0.01<LHR<0.1
    # row 15 -- ONI=warm, LKB<=1Mt, LHR>=0.1
    # row 16 -- ONI=warm, LKB>1Mt, LHR<=0.01
    # row 17 -- ONI=warm, LKB<=1Mt, 0.01<LHR<0.1
    # row 18 -- ONI=warm, LKB>1Mt, LHR>=0.1
    for(i in 1:18){
      mu.index.new[i]<-inprod(predX[i,],beta[]) # posterior expectation at new data points
      index.new[i]~dnorm(mu.index.new[i],tau.index) # posterior predictive
    }
    
    # some interesting probabilities
    
    # that effects change expected performance relative to the reference case
    # high biomass
    prob[1]<-ifelse(mu.index.new[2]<mu.index.new[1],1,0)
    prob.new[1]<-ifelse(index.new[2]<index.new[1],1,0)
    # med hr
    prob[2]<-ifelse(mu.index.new[3]<mu.index.new[1],1,0)
    prob.new[2]<-ifelse(index.new[3]<index.new[1],1,0)
    # high hr
    prob[3]<-ifelse(mu.index.new[5]<mu.index.new[1],1,0)
    prob.new[3]<-ifelse(index.new[5]<index.new[1],1,0)
    # neutral ONI
    prob[4]<-ifelse(mu.index.new[7]<mu.index.new[1],1,0)
    prob.new[4]<-ifelse(index.new[7]<index.new[1],1,0)
    # warm ONI
    prob[5]<-ifelse(mu.index.new[13]<mu.index.new[1],1,0)
    prob.new[5]<-ifelse(index.new[13]<index.new[1],1,0)
    # worst case
    prob[6]<-ifelse(mu.index.new[12]<mu.index.new[1],1,0)
    prob.new[6]<-ifelse(index.new[12]<index.new[1],1,0)
    
    # that other effects are more extreme than environmental effects
    # med hr has more negative effect than neutral ONI
    prob[7]<-ifelse(mu.index.new[3]<mu.index.new[7],1,0)
    prob.new[7]<-ifelse(index.new[3]<index.new[7],1,0)
    # that high hr has more negative effect than neutral ONI
    prob[8]<-ifelse(mu.index.new[5]<mu.index.new[7],1,0)
    prob.new[8]<-ifelse(index.new[5]<index.new[7],1,0)
    # that high krill biomass has more negative effect than neutral ONI
    prob[9]<-ifelse(mu.index.new[2]<mu.index.new[7],1,0)
    prob.new[9]<-ifelse(index.new[2]<index.new[7],1,0)
    # that med hr has more negative effect than warm ONI
    prob[10]<-ifelse(mu.index.new[3]<mu.index.new[13],1,0)
    prob.new[10]<-ifelse(index.new[3]<index.new[13],1,0)
    # that high hr has more negative effect than warm ONI
    prob[11]<-ifelse(mu.index.new[5]<mu.index.new[13],1,0)
    prob.new[11]<-ifelse(index.new[5]<index.new[13],1,0)
    # that high krill biomass has more negative effect than warm ONI
    prob[12]<-ifelse(mu.index.new[2]<mu.index.new[13],1,0)
    prob.new[12]<-ifelse(index.new[2]<index.new[13],1,0)
    
    
    # that effects change expected performance relative to long-term mean
    # reference case
    prob[13]<-ifelse(mu.index.new[1]<0,1,0)
    prob.new[13]<-ifelse(index.new[1]<0,1,0)
    # high biomass
    prob[14]<-ifelse(mu.index.new[2]<0,1,0)
    prob.new[14]<-ifelse(index.new[2]<0,1,0)
    # med hr
    prob[15]<-ifelse(mu.index.new[3]<0,1,0)
    prob.new[15]<-ifelse(index.new[3]<0,1,0)
    # high hr
    prob[16]<-ifelse(mu.index.new[5]<0,1,0)
    prob.new[16]<-ifelse(index.new[5]<0,1,0)
    # neutral ONI
    prob[17]<-ifelse(mu.index.new[7]<0,1,0)
    prob.new[17]<-ifelse(index.new[7]<0,1,0)
    # warm ONI
    prob[18]<-ifelse(mu.index.new[13]<0,1,0)
    prob.new[18]<-ifelse(index.new[13]<0,1,0)
    # worst case
    prob[19]<-ifelse(mu.index.new[12]<0,1,0)
    prob.new[19]<-ifelse(index.new[12]<0,1,0)

  }
"


# objects needed to fit the model and monitor variables of interest
# there's a trick here -- if is.na(survey) then make survey a big number to prevent
# division by zero during imputation procedure these will either be replaced
# by imputed values (summer surveys) or not used (winter surveys)
#

pred.matrix<-matrix(c(1,-1,-1,-1,-1,-1,
                      1,1,-1,-1,-1,-1,
                      1,-1,1,0,-1,-1,
                      1,1,1,0,-1,-1,
                      1,-1,0,1,-1,-1,
                      1,1,0,1,-1,-1,
                      
                      1,-1,-1,-1,1,0,
                      1,1,-1,-1,1,0,
                      1,-1,1,0,1,0,
                      1,1,1,0,1,0,
                      1,-1,0,1,1,0,
                      1,1,0,1,1,0,
                      
                      1,-1,-1,-1,0,1,
                      1,1,-1,-1,0,1,
                      1,-1,1,0,0,1,
                      1,1,1,0,0,1,
                      1,-1,0,1,0,1,
                      1,1,0,1,0,1),nrow=18,ncol=6,byrow=TRUE)

hr.data<-list(index=as.vector(junk$index),
              bclass=junk$bclass,
              hrclass=junk$hrclass,
              oniclass=as.numeric(factor(junk$oni.class)),
              nobs=dim(junk)[1],
              predX=pred.matrix)

hr.params<-c("beta","sd.index","t.sd.index")

beta.init1<-rep(-1,6)
beta.init2<-rep(0,6)
beta.init3<-rep(1,6)


hr.inits<-list(list(beta=beta.init1,t.sd.index=0.1,.RNG.seed=123,.RNG.name="base::Super-Duper"),
               list(beta=beta.init2,t.sd.index=1.0,.RNG.seed=456,.RNG.name="base::Super-Duper"),
               list(beta=beta.init3,t.sd.index=1.9,.RNG.seed=789,.RNG.name="base::Super-Duper"))

hr.derived<-c("index.new","mu.index.new","prob","prob.new")


# now do the analysis
library(coda)
library(rjags)

hr.jags<-jags.model(textConnection(modelstring),hr.data,hr.inits,n.chains=3,n.adapt=50000)
# burn in for 150000 iterations
update(hr.jags, n.iter=100000)
hr.params.post<-coda.samples(hr.jags,hr.params,n.iter=50000,thin=25)
hr.derived.post<-coda.samples(hr.jags,hr.derived,n.iter=50000,thin=25)
hr.params.summ<-summary(hr.params.post)
hr.derived.summ<-summary(hr.derived.post)

# write input/output
out.dir <- "./A_Results/Model02_SSMU/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(hr.params.post, file.path(out.dir, file = "hr.params.post.rds"))
saveRDS(hr.derived.post, file.path(out.dir, file = "hr.derived.post.rds"))

saveRDS(hr.params.summ, file.path(out.dir, file = "hr.params.summ.rds"))
saveRDS(hr.derived.summ, file.path(out.dir, file = "hr.derived.summ.rds"))


### Step 3: PLOTTING  ###

require(ggmcmc)
hr.params.s<-ggs(hr.params.post)
hr.derived.s<-ggs(hr.derived.post)
# 
# # just want to copy hr.params.s to work with it for plotting diagnostics without screwing up the original object
# # also get rid of chains for t.sd.index since this is not really a parameter of interest
# HR.labels<-data.frame(Parameter=dimnames(hr.params.post[[1]])[[2]],
#                       Label=c("alpha","beta[3]","beta[4]","beta[5]","beta[1]","beta[2]",
#                               "K[B,-]","K[D,-]","K[B,+]","K[D,+]","sigma","phi","exclude"))

# Modified, correct HR.labels:
HR.labels <- data.frame(
  Parameter = dimnames(hr.params.post[[1]])[[2]],
  Label = c("alpha", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", 
            "sigma", "phi")
) 
hr.params2.s<-ggs(hr.params.post,par_labels = HR.labels)
hr.params2.s<-hr.params2.s[hr.params2.s$ParameterOriginal!="t.sd.index",]




# Revised Supplementary Figure S8 (Watters et al.)

# reference (best case)
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-1.75,1.75),xaxt="n",xlim=c(0.5,7.5),
        ylab="expected performance", xlab="",whisklty=1,boxwex=1,at=1,col="white")
# ONI
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(25,31))),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=2:3,col="gray80")
# biomass
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==20),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=1,add=TRUE,at=4,col="gray40",medcol="white")
# harvest rate
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(21,23))),
        range=0,yaxt="n",xaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=5:6,col="black",medcol="white")
# worst case
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==36),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=1,add=TRUE,at=7,col="white")
axis(1,at=1:7,labels=c("best case","-0.5 < ONI < 0.5","ONI >= 0.5","LKB > 1 Mt",
                       "0.01< LHR <0.10","LHR >= 0.10","worst case"), cex=0.75)
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)
## ADD BLUE LINE
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added

dev.copy(
  png,
  filename = file.path(out.dir, file = "Penguin_performance_mod02_SSMU.png"),
  width    = 12,
  height   = 7,
  units    = "in",
  res      = 300
)
dev.off()

write.csv(hr.derived.s, file = file.path(out.dir, file ="hr_derived_s.csv"))

# # ***********************************************************************************************
# # *************************************** ADDED FIGURE - plot all cases***************************
# # ***********************************************************************************************

# # reference (best case)
# Set up plot margins to accommodate rotated labels (bottom margin increased)
par(mar=c(8.5, 4, 4, 0.5) + 0.1)
# Color ONI classes
case_colors <- c(rep("lightblue", 6), rep("white", 6), rep("lightpink", 6))
# Create x-axis labels
x_labels <- c("LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1",
              "LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1",
              "LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1")

boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,
        subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,18.5),
        ylab="expected performance",xlab="",
        whisklty=1,boxwex=1,at=1,col=case_colors[1])

for(i in 2:18){
  boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,
          subset=(as.numeric(hr.derived.s$Parameter)==(i+18)),
          range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,18.5),
          whisklty=1,boxwex=1,at=i,
        col=case_colors[i],add=TRUE)
}

# axis(1,at=1:18,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))
axis(1, at=1:18, labels=FALSE)
text(x=1:18, y=par("usr")[3]-0.3, labels=x_labels, srt=90, adj=1, xpd=TRUE, cex=0.65)

abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added

# Add legend for ONI classes
legend("topright", 
       legend=c("ONI < -0.5C", "-0.5C < ONI < 0.5C", "ONI > 0.5C"),
       fill=c("lightblue", "white", "lightpink"),
       border="black",
       title="ONI Classes",
       cex=0.65)

dev.copy(
  png,
  filename = file.path(out.dir, file = "Penguin_performance_ALL_mod02_SSMU.png"),
  width    = 8,
  height   = 6,
  units    = "in",
  res      = 300
)
dev.off()

# ****************************************************************************************************

# Supplementary Figures

# ****************************************************************************************************

# S1 -- trace plots of main model parameters
ggs_traceplot(hr.params2.s) + facet_wrap(~ Parameter, ncol = 3, scales="free")

# S2 -- scale-reduction factors
ggs_Rhat(hr.params2.s)

# S3 -- Geweke Z-scores
ggs_geweke(hr.params2.s,shadow_limit = 1.96)

# S4 -- autocorrelation plots
ggs_autocorrelation(hr.params2.s)

# S5 -- crosscorrelations
ggs_crosscorrelation(hr.params2.s)
#hr.params2.s<-ggs(hr.params.post)
#hr.params2.s<-hr.params2.s[hr.params2.s$Parameter!="t.sd.index",]
#ggs_pairs(hr.params2.s,lower=list(continuous="density"))

# S6 -- posterior distributions
ggs_histogram(hr.params2.s) + facet_wrap(~ Parameter, ncol = 3, scales="free")
#ggs_density(hr.params2.s) + facet_wrap(~ Parameter, ncol = 3, scales="free")

# S7 -- plot posterior predictive distributions over data for visual posterior predictive check
xx<-junk
xx$case<-ifelse(xx$oniclass==1 & xx$bclass==1 & xx$hrclass==1,1,
          ifelse(xx$oniclass==1 & xx$bclass==2 & xx$hrclass==1,2,
           ifelse(xx$oniclass==1 & xx$bclass==1 & xx$hrclass==2,3,
             ifelse(xx$oniclass==1 & xx$bclass==2 & xx$hrclass==2,4,
              ifelse(xx$oniclass==1 & xx$bclass==1 & xx$hrclass==3,5,
                ifelse(xx$oniclass==1 & xx$bclass==2 & xx$hrclass==3,6,
                ifelse(xx$oniclass==2 & xx$bclass==1 & xx$hrclass==1,7,
                  ifelse(xx$oniclass==2 & xx$bclass==2 & xx$hrclass==1,8,
                    ifelse(xx$oniclass==2 & xx$bclass==1 & xx$hrclass==2,9,
                      ifelse(xx$oniclass==2 & xx$bclass==2 & xx$hrclass==2,10,
                        ifelse(xx$oniclass==2 & xx$bclass==1 & xx$hrclass==3,11,
                          ifelse(xx$oniclass==2 & xx$bclass==2 & xx$hrclass==3,12,
                            ifelse(xx$oniclass==3 & xx$bclass==1 & xx$hrclass==1,13,
                              ifelse(xx$oniclass==3 & xx$bclass==2 & xx$hrclass==1,14,
                                ifelse(xx$oniclass==3 & xx$bclass==1 & xx$hrclass==2,15,
                                  ifelse(xx$oniclass==3 & xx$bclass==2 & xx$hrclass==2,16,
                                    ifelse(xx$oniclass==3 & xx$bclass==1 & xx$hrclass==3,17,18)))))))))))))))))

#### Save 'xx'
saveRDS(xx, file.path(out.dir, "Data_Fig_S7.rds"))  # save the data
####

# Create mapping for colors and hatching based on case number
# Cases 1-6: oniclass=1 (light blue)
# Cases 7-12: oniclass=2 (transparent/white)
# Cases 13-18: oniclass=3 (light red)
case_colors <- c(rep("lightblue", 6), rep("white", 6), rep("lightpink", 6))

# Create x-axis labels
x_labels <- c("LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1",
              "LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1",
              "LKB<=1Mt, LHR<=0.01", "LKB>1Mt, LHR<=0.01",
              "LKB<=1Mt, 0.01<LHR<0.1", "LKB>1Mt, 0.01<LHR<0.1",
              "LKB<=1Mt, LHR>=0.1", "LKB>1Mt, LHR>=0.1")

# Set up plot margins to accommodate rotated labels (bottom margin increased)
par(mar=c(8.5, 4, 4, 0.5) + 0.1)

# plot the data
plot(jitter(xx$case,amount=0.25),xx$index,type="n",xlim=c(0.65,18.35),ylim=c(-3.5,3.5),
     xlab="",xaxt="n",ylab="Performance index",pch=16)
# add posterior predictive distributions as box plots
for(i in 1:18){
  boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,
              subset=(as.numeric(hr.derived.s$Parameter)==i),
              range=1.5,outline=FALSE,xaxt="n",whisklty=1,boxwex=1,at=i,
              col=case_colors[i],add=TRUE)
}

# plot the data
points(jitter(xx$case[xx$season=="S"],amount=0.2),xx$index[xx$season=="S"],cex=0.5,col="red",pch=16)
points(jitter(xx$case[xx$season=="W"],amount=0.2),xx$index[xx$season=="W"],cex=0.5,col="blue",pch=16)
abline(h=0)
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
## ADD BLUE LINE
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added

# Add custom x-axis with rotated labels
axis(1, at=1:18, labels=FALSE)
text(x=1:18, y=par("usr")[3]-0.3, labels=x_labels, srt=90, adj=1, xpd=TRUE, cex=0.65)

# Add legend for ONI classes
legend("topright", 
       legend=c("ONI < -0.5C", "-0.5C < ONI < 0.5C", "ONI > 0.5C"),
       fill=c("lightblue", "white", "lightpink"),
       border="black",
       title="ONI Classes",
       cex=0.65)

dev.copy(
  png,
  filename = file.path(out.dir, file = "Posterior_predictive_dist_mod02_SSMU.png"),
  width    = 8,
  height   = 6,
  units    = "in",
  res      = 300
)
dev.off()

##########################################################################

# ============================================================
# Observed data plots — Mod_02_SSMUs
# ============================================================

# Map individual SSMUs back to gSSMU for faceting (consistent with previous plots)
# gSSMU 1 (Bransfield): APBSE, APBSW
# gSSMU 2 (Drake Passage): APDPE, APDPW, APEI

out.dir.mod02 <- "./A_Results/Model02_SSMU/"

junk_plot <- junk %>%
  mutate(
    gSSMU_lab = case_when(
      SSMU %in% c("APBSE", "APBSW")       ~ "gSSMU 1 (Bransfield)",
      SSMU %in% c("APDPE", "APDPW", "APEI") ~ "gSSMU 2 (Drake Passage)",
      TRUE ~ "Other"
    ),
    sam_lab     = ifelse(sam.sign == "Neg", "SAM Negative", "SAM Positive"),
    hr          = catch / survey,
    bclass_lab  = factor(ifelse(bclass == 1, "≤ 1 Mt", "> 1 Mt"),
                         levels = c("≤ 1 Mt", "> 1 Mt")),
    hrclass_lab = factor(case_when(
      hrclass == 1 ~ "≤ 0.01",
      hrclass == 2 ~ "0.01 – <0.10",
      hrclass == 3 ~ "≥ 0.10"),
      levels = c("≤ 0.01", "0.01 – <0.10", "≥ 0.10"))
  )

# Summer-only subset (for p1, p1b, p2)
summer_obs02  <- junk_plot %>% filter(season == "S")
n_summer_obs02 <- nrow(summer_obs02)
n_all_obs02    <- nrow(junk_plot)

# ---- p1: summer biomass ----
p1_obs02 <- ggplot(summer_obs02, aes(x = survey / 1e6,
                                     y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "#66bb6a", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 0.15)) +
  scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  labs(
    title    = "Observed summer krill biomass — Model SSMUs",
    x        = "Biomass (million tonnes)",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p1_obs02

# ---- p1b: log biomass ----
p1b_obs02 <- ggplot(summer_obs02, aes(x = log(survey),
                                      y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "#2d7a3d", color = "white", alpha = 0.85) +
  geom_vline(xintercept = log(1e6), linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(
    title    = "Observed summer krill biomass (log scale) — Model SSMUs",
    x        = "ln(Biomass in tonnes)",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p1b_obs02

# ---- p2: summer harvest rate ----
p2_obs02 <- ggplot(summer_obs02, aes(x = hr,
                                     y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 0.005, boundary = 0, fill = "#66bb6a", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 0.10, linetype = "solid",  color = "black", linewidth = 0.7) +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.30)) +
  scale_x_continuous(breaks = seq(0, 0.25, by = 0.01)) +
  labs(
    title    = "Observed summer harvest rate — Model SSMUs",
    x        = "Harvest Rate",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p2_obs02

# ---- p3a: biomass class (all observations) ----
p3a_obs02 <- summer_obs02 %>%
  count(gSSMU_lab, sam_lab, bclass_lab) %>%
  ggplot(aes(x = bclass_lab, y = n/sum(n), fill = bclass_lab)) +
  geom_col(width = 0.55, color = "white") +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(ylim = c(0, 0.30)) +
  scale_fill_manual(values = c("≤ 1 Mt" = "#a8d5ba", "> 1 Mt" = "#2d7a3d")) +
  labs(
    title    = "Biomass class distribution — Model SSMUs",
    x        = "Biomass class",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none"); p3a_obs02

# ---- p3b: harvest rate class (all observations) ----
p3b_obs02 <- summer_obs02 %>%
  count(gSSMU_lab, sam_lab, hrclass_lab) %>%
  ggplot(aes(x = hrclass_lab, y = n/sum(n), fill = hrclass_lab)) +
  geom_col(width = 0.55, color = "white") +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(ylim = c(0, 0.30)) +
  scale_fill_manual(values = c(
    "≤ 0.01"       = "#a8d5ba",
    "0.01 – <0.10" = "#66bb6a",
    "≥ 0.10"       = "#2d7a3d"
  )) +
  labs(
    title    = "Harvest rate class distribution — Model SSMUs",
    x        = "Harvest rate class",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none"); p3b_obs02


p3b_obs02 <- summer_obs02 %>%
  mutate(hrclass_lab = factor(hrclass_lab, 
                              levels = c("≤ 0.01", "0.01 – <0.10", "≥ 0.10"))) %>%
  count(gSSMU_lab, sam_lab, hrclass_lab, .drop = FALSE) %>%
  ggplot(aes(x = hrclass_lab, y = n / sum(n), fill = hrclass_lab)) +
  geom_col(width = 0.55, color = "white") +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(ylim = c(0, 0.30)) +
  scale_fill_manual(values = c(
    "≤ 0.01"       = "#a8d5ba",
    "0.01 – <0.10" = "#66bb6a",
    "≥ 0.10"       = "#2d7a3d"
  )) +
  labs(
    title    = "Harvest rate class distribution — Model SSMUs",
    x        = "Harvest rate class",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none"); p3b_obs02


# ---- save ----
ggsave(file.path(out.dir.mod02, "obs_freq_dist_summer_biomass.png"),     p1_obs02,  width=8, height=6, dpi=300)
ggsave(file.path(out.dir.mod02, "obs_freq_dist_summer_biomass_log.png"), p1b_obs02, width=8, height=6, dpi=300)
ggsave(file.path(out.dir.mod02, "obs_freq_dist_summer_hr.png"),          p2_obs02,  width=8, height=6, dpi=300)
ggsave(file.path(out.dir.mod02, "obs_freq_dist_bclass.png"),             p3a_obs02, width=8, height=6, dpi=300)
ggsave(file.path(out.dir.mod02, "obs_freq_dist_hrclass.png"),            p3b_obs02, width=8, height=6, dpi=300)


