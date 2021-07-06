library(tidyverse)
dir = "E:/Documents/R workspace/Southern Ocean/Watters Sci Rep/Data from paper"
setwd(dir)
make.localhr.data<-function(trim=1,plot.winter=FALSE){
  ###########################################################################################################
  # generate the summer indices
  ###########################################################################################################
  #
  # fledge weight (fwt)
  # bigger indicates better summer
  fwt<-read.csv("fweight.csv",header=TRUE,stringsAsFactors = FALSE)
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
  phs<-read.csv("success.csv",header=TRUE,stringsAsFactors = FALSE)
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
  td<-read.csv("tripduration.csv",header=TRUE,stringsAsFactors = FALSE)
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
  ade1<-read.csv("massatlay.csv",header=TRUE,stringsAsFactors = FALSE)
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
  e1e2<-read.csv("egg.csv",header=TRUE,stringsAsFactors = FALSE)
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
  cid<-read.csv("cid.csv",header=TRUE,stringsAsFactors = FALSE)[,1:4]
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
  rec<-read.csv("recruitment.csv",header=TRUE,stringsAsFactors = FALSE)[,1:4]
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
  survey<-read.csv("krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)
  # use next line if want to filter acoustic data to have minimum number of miles (comment out if not desired)
  # as per CSR, 80 nmi would be about equivalent of 2 tracklines in the Bransfield
  #survey<-survey[survey$nmi.count>=80,]
  # could try changing "biomass" in following line to "mean.density.gm2" or "median.density.gm2" but haven't done that
  # change here
  survey.mod = as_tibble(survey) %>%
    group_split(gSSMU)
  # 28.7/(28.7+22) #fraction of nmi area # change here
  survey.mod[[1]]$biomass = survey.mod[[1]]$biomass * 0.566075 
  # 15.8/(15.8+16.4+36.2)fraction of nmi area # change here 
  survey.mod[[2]]$biomass = survey.mod[[2]]$biomass * 0.230994 
  # 22.0/(28.7+22.0) change here
  survey.mod[[3]]$biomass = survey.mod[[3]]$biomass * 0.433925 
  # change here
  survey.mod[[4]]$biomass = survey.mod[[4]]$biomass * 1 
  # change here
  survey = data.frame(bind_rows(survey.mod[[1]],survey.mod[[2]],survey.mod[[3]],survey.mod[[4]]))  
  # use next line if want to filter acoustic data to have minimum number of miles (comment out if not desired)
  # as per CSR, 80 nmi would be about equivalent of 2 tracklines in the Bransfield
  #survey<-survey[survey$nmi.count>=80,]
  # could try changing "biomass" in following line to "mean.density.gm2" or "median.density.gm2" but haven't done that
  survey<-tapply(survey$biomass,list(survey$Year,survey$gSSMU),mean,na.rm=TRUE)
  survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                     gSSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                     survey=c(survey),stringsAsFactors = FALSE)
  survey$season<-ifelse(survey$cal.yr<2012,"S","W")
  # use next line if want to remove winter survey data altogether (comment out if not desired)
  #survey<-survey[survey$season=="S",]
  survey$matchme<-paste(survey$cal.yr,survey$season,survey$gSSMU,sep="|")
  #print(str(survey))
  #
  ###########################################################################################################
  #
  # krill fishery catches
  fishery<-read.csv("c1.csv",header=TRUE,stringsAsFactors = FALSE)
  fishery$season<-ifelse(is.element(fishery$Month,c(10:12,1:2)),"S","W") #1:3 if March in Summer
  gSSMU1<-c("APBSE") # change here
  gSSMU2<-c("APDPW") # change here
  gSSMU3<-"APBSW" # change here
  gSSMU4<-c("APW","APE") # change here
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
  out$gSSMU<-rep(NA,dim(out)[1])
  out$gSSMU<-ifelse(out$SPECIES=="ADPE"&out$PROJECT=="COPA"&out$season=="S",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="ADPE"&out$PROJECT=="COPA"&out$season=="W",NA,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="S",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="W",NA,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="S",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="W",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="S",2,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="W",NA,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="S",2,out$gSSMU)
  # use following line if GEPE at CS forage in gSSMU 2 during winter
  #out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="W",2,out$gSSMU)
  # use following line if GEPE at CS forage in gSSMU 1 during winter
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="CS"&out$season=="W",1,out$gSSMU)
  #
  out$matchme<-paste(out$cal.yr,out$season,out$gSSMU,sep="|")
  out$survey<-survey$survey[match(out$matchme,survey$matchme)]
  out$catch<-fishery$catch[match(out$matchme,fishery$matchme)]
  #
  out<-out[!is.na(out$gSSMU),]
  #
  ###########################################################################################################
  # pull in the environmental indices
  ###########################################################################################################
  #
  # SOUTHERN ANNULAR MODE
  sam<-read.csv("sam.csv")
  names(sam)<-c("yr","mo","sam")
  sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
  sam$YEAR<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
  sam<-tapply(sam$sam,list(sam$YEAR,sam$season),mean)
  sam<-data.frame(YEAR=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))
  out$sam<-sam$sam[match(paste(out$cal.yr,out$season,sep="|"),paste(sam$YEAR,sam$season,sep="|"))]
  out$sam.sign<-ifelse(out$sam<0,"Neg","Pos")
  #
  # OCEANIC NINO INDEX
  oni<-read.csv("oni.csv",stringsAsFactors = FALSE)
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

junk<-junk[!is.na(junk$survey),]

junk$bclass<-ifelse(junk$survey<=1000000,1,2)
junk$hrclass<-ifelse(junk$catch/junk$survey<=0.01,1,ifelse(junk$catch/junk$survey>=0.1,3,2))
junk$oniclass<-ifelse(junk$oni.class=="Cool",1,ifelse(junk$oni.class=="Warm",3,2))




##################################################################################################




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
              oniclass=junk$oniclass,
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
hr.params.post<-coda.samples(hr.jags,hr.params,n.iter=125000,thin=25)
hr.derived.post<-coda.samples(hr.jags,hr.derived,n.iter=125000,thin=25)
hr.params.summ<-summary(hr.params.post)
hr.derived.summ<-summary(hr.derived.post)


require(ggmcmc)
hr.params.s<-ggs(hr.params.post)
hr.derived.s<-ggs(hr.derived.post)

save(hr.derived.s, file="E:/Documents/R workspace/Southern Ocean/Watters Sci Rep/Watters and Kruger/WattKrug/hr_derived_s_37_nonimpute.csv")
#ggmcmc(hr.params.s,file="diagnostics_hr.pdf")

# ************************************************************************************

# plot posterior expectations of marginal effects

# ************************************************************************************


# Supplementary Figure S8
# reference
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-1.75,1.75),xaxt="n",xlim=c(0.5,7.5),ylab="expected performance",whisklty=1,boxwex=1,at=1)
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
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=1,add=TRUE,at=7)
axis(1,at=1:7,labels=c("reference","-0.5 < ONI < 0.5","ONI >= 0.5","LKB > 1 Mt","0.01 < LHR < 0.10","LHR >= 0.10","worst case"))
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)





# ****************************************************************************************************

# plot posterior predictive distributions overlapping some data for a quick posterior predictive check

# ****************************************************************************************************


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
# plot the data
plot(jitter(xx$case,amount=0.25),xx$index,type="n",xlim=c(0.65,18.35),ylim=c(-3.5,3.5),xlab="case",xaxt="n",ylab="std performance index",pch=16)
# add posterior predictive distributions as box plots
for(i in 1:18){
  boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==i),
          range=1.5,outline=FALSE,xaxt="n",whisklty=1,boxwex=1,at=i,add=TRUE)
}
# plot the data
points(jitter(xx$case[xx$season=="S"],amount=0.2),xx$index[xx$season=="S"],cex=0.5,col="red",pch=16)
points(jitter(xx$case[xx$season=="W"],amount=0.2),xx$index[xx$season=="W"],cex=0.5,col="blue",pch=16)
axis(1,at=1:18)
abline(h=0)
