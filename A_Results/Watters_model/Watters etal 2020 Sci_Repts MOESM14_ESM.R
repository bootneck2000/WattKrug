# Original Code - Watters et al (2020)
# summer LKB IMPUTED

library(tidyverse)
library(lattice)

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
  out$gSSMU<-ifelse(out$SPECIES=="ADPE"&out$PROJECT=="COPA"&out$season=="W",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="S",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="COPA"&out$season=="W",2,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="S",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="GEPE"&out$PROJECT=="COPA"&out$season=="W",1,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="S",2,out$gSSMU)
  out$gSSMU<-ifelse(out$SPECIES=="CHPE"&out$PROJECT=="CS"&out$season=="W",2,out$gSSMU)
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


##################################################################################################


modelstring<-"
  # George Watters -- April 2019
  
  # missing estimates of krill biomass during the summer are imputed based on observed relationship
  # with environmental covariates and the uncertainty in this imputation is carried forward into
  # estimation of effects from krill availability

  # samclass are environment classes based on SAM where
  # 1 = negative
  # 2 = positive
    
  # oniclass are environment classes based on ONI where
  # 1 = cool
  # 2 = neutral
  # 3 = warm

  # studies are combinations of site, species, and monitoring parameter
    
  

  model{

    #####################################################################
    #
    # model for missing biomass estimates in summer (do not think data are sufficient to do this for winter)
    # this idea is based on example at http://www.columbia.edu/~cjd11/charles_dimaggio/DIRE/styled-4/styled-11/code-10/
    #
    #####################################################################    

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

   

    #####################################################################
    #
    # model for penguin performance
    #
    #####################################################################


    # create the design matrix with sum-to-zero constraints
    for(i in 1:nobs){
      X[i,1]<-1.0     # intercept
      X[i,2]<-equals(bclass[i],2)-equals(bclass[i],1) # b2
      X[i,3]<-equals(hrclass[i],2)-equals(hrclass[i],1) # hr2
      X[i,4]<-equals(hrclass[i],3)-equals(hrclass[i],1) # hr3
      X[i,5]<-equals(oniclass[i],2)-equals(oniclass[i],1) # o2
      X[i,6]<-equals(oniclass[i],3)-equals(oniclass[i],1) # o3
    }
  

    # first the likelihood
    # loop over number of data points (nobs)
    for(i in 1:nobs){
      index[i]~dnorm(mu[i],tau.index)
      mu[i] <- inprod(X[i,],beta[])
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
    # row 1 -- ONI=cool, LKB<=1Mt, LHR<=0.01 (reference or best case)
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
      mu.new[i]<-inprod(predX[i,],beta[]) # posterior expectation at new data points
      index.new[i]~dnorm(mu.new[i],tau.index) # posterior predictive
    }
    
    # some interesting probabilities
    
    # that effects change expected performance relative to the reference case
    # high biomass
    prob[1]<-ifelse(mu.new[2]<mu.new[1],1,0)
    prob.new[1]<-ifelse(index.new[2]<index.new[1],1,0)
    # med hr
    prob[2]<-ifelse(mu.new[3]<mu.new[1],1,0)
    prob.new[2]<-ifelse(index.new[3]<index.new[1],1,0)
    # high hr
    prob[3]<-ifelse(mu.new[5]<mu.new[1],1,0)
    prob.new[3]<-ifelse(index.new[5]<index.new[1],1,0)
    # neutral ONI
    prob[4]<-ifelse(mu.new[7]<mu.new[1],1,0)
    prob.new[4]<-ifelse(index.new[7]<index.new[1],1,0)
    # warm ONI
    prob[5]<-ifelse(mu.new[13]<mu.new[1],1,0)
    prob.new[5]<-ifelse(index.new[13]<index.new[1],1,0)
    # worst case
    prob[6]<-ifelse(mu.new[12]<mu.new[1],1,0)
    prob.new[6]<-ifelse(index.new[12]<index.new[1],1,0)
    
    # that other effects are more extreme than environmental effects
    # med hr has more negative effect than neutral ONI
    prob[7]<-ifelse(mu.new[3]<mu.new[7],1,0)
    prob.new[7]<-ifelse(index.new[3]<index.new[7],1,0)
    # that high hr has more negative effect than neutral ONI
    prob[8]<-ifelse(mu.new[5]<mu.new[7],1,0)
    prob.new[8]<-ifelse(index.new[5]<index.new[7],1,0)
    # that high krill biomass has more negative effect than neutral ONI
    prob[9]<-ifelse(mu.new[2]<mu.new[7],1,0)
    prob.new[9]<-ifelse(index.new[2]<index.new[7],1,0)
    # that med hr has more negative effect than warm ONI
    prob[10]<-ifelse(mu.new[3]<mu.new[13],1,0)
    prob.new[10]<-ifelse(index.new[3]<index.new[13],1,0)
    # that high hr has more negative effect than warm ONI
    prob[11]<-ifelse(mu.new[5]<mu.new[13],1,0)
    prob.new[11]<-ifelse(index.new[5]<index.new[13],1,0)
    # that high krill biomass has more negative effect than warm ONI
    prob[12]<-ifelse(mu.new[2]<mu.new[13],1,0)
    prob.new[12]<-ifelse(index.new[2]<index.new[13],1,0)
    
    
    # that effects change expected performance relative to long-term mean
    # reference case
    prob[13]<-ifelse(mu.new[1]<0,1,0)
    prob.new[13]<-ifelse(index.new[1]<0,1,0)
    # high biomass
    prob[14]<-ifelse(mu.new[2]<0,1,0)
    prob.new[14]<-ifelse(index.new[2]<0,1,0)
    # med hr
    prob[15]<-ifelse(mu.new[3]<0,1,0)
    prob.new[15]<-ifelse(index.new[3]<0,1,0)
    # high hr
    prob[16]<-ifelse(mu.new[5]<0,1,0)
    prob.new[16]<-ifelse(index.new[5]<0,1,0)
    # neutral ONI
    prob[17]<-ifelse(mu.new[7]<0,1,0)
    prob.new[17]<-ifelse(index.new[7]<0,1,0)
    # warm ONI
    prob[18]<-ifelse(mu.new[13]<0,1,0)
    prob.new[18]<-ifelse(index.new[13]<0,1,0)
    # worst case
    prob[19]<-ifelse(mu.new[12]<0,1,0)
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
              survey=ifelse(is.na(junk$survey),1E12,junk$survey),
              catch=junk$catch,
              gssmu=junk$gSSMU,
              oniclass=as.numeric(factor(junk$oni.class)),
              samclass=as.numeric(factor(junk$sam.sign)),
              summer=junk$survey[junk$season=="S"],
              impute.me=ifelse(is.na(junk$survey)&junk$season=="S",1,0),
              meanlogsummer=tapply(log(junk$survey[junk$season=="S"]),list(junk$gSSMU[junk$season=="S"],junk$sam.sign[junk$season=="S"]),mean,na.rm=TRUE),
              sdlogsummer=sd(log(junk$survey[junk$season=="S"]),na.rm=TRUE),
              nobs=dim(junk)[1],
              nsummerobs=as.vector(table(junk$season)[1]),
              predX=pred.matrix)

hr.params<-c("beta","mulogsummer","sigmalogsummer","sd.index","t.sd.index")

beta.init1<-rep(-1,6)
beta.init2<-rep(0,6)
beta.init3<-rep(1,6)


hr.inits<-list(list(beta=beta.init1,t.sd.index=0.1,.RNG.seed=123,.RNG.name="base::Super-Duper"),
                list(beta=beta.init2,t.sd.index=1.0,.RNG.seed=456,.RNG.name="base::Super-Duper"),
                list(beta=beta.init3,t.sd.index=1.9,.RNG.seed=789,.RNG.name="base::Super-Duper"))

hr.derived<-c("index.new","mu.new","prob","prob.new")

hr.imputed<-"hr"



# now do the analysis
library(coda)
library(rjags)

hr.jags<-jags.model(textConnection(modelstring),hr.data,hr.inits,n.chains=3,n.adapt=50000)
# burn in for 150000 iterations
update(hr.jags, n.iter=100000)
hr.params.post<-coda.samples(hr.jags,hr.params,n.iter=50000,thin=25)
hr.derived.post<-coda.samples(hr.jags,hr.derived,n.iter=50000,thin=25)
hr.imputed.post<-coda.samples(hr.jags,hr.imputed,n.iter=50000,thin=25)
hr.params.summ<-summary(hr.params.post)
hr.derived.summ<-summary(hr.derived.post)
hr.imputed.summ<-summary(hr.imputed.post)

require(ggmcmc)
hr.params.s<-ggs(hr.params.post)
hr.derived.s<-ggs(hr.derived.post)

# just want to copy hr.params.s to work with it for plotting diagnostics without screwing up the original object
# also get rid of chains for t.sd.index since this is not really a parameter of interest
HR.labels<-data.frame(Parameter=dimnames(hr.params.post[[1]])[[2]],
                      Label=c("alpha","beta[3]","beta[4]","beta[5]","beta[1]","beta[2]","K[B,-]","K[D,-]","K[B,+]","K[D,+]","sigma","phi","exclude"))

hr.params2.s<-ggs(hr.params.post,par_labels = HR.labels)
hr.params2.s<-hr.params2.s[hr.params2.s$ParameterOriginal!="t.sd.index",]

#ggmcmc(hr.params.s,file="diagnostics_hr_params_final.pdf")
#ggmcmc(hr.derived.s,file="diagnostics_hr_derived_final.pdf")

# write input/output
out.dir <- "./A_Results/Watters_model/Imputed/"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(hr.params.post, file.path(out.dir, file = "hr.params.post.rds"))
saveRDS(hr.derived.post, file.path(out.dir, file = "hr.derived.post.rds"))

saveRDS(hr.params.summ, file.path(out.dir, file = "hr.params.summ.rds"))
saveRDS(hr.derived.summ, file.path(out.dir, file = "hr.derived.summ.rds"))


# ************************************************************************************

# plot posterior expectations of marginal effects

# ************************************************************************************

# Figure 2

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
# axis(1,at=1:7,labels=c("best case","-0.5 < ONI < 0.5","ONI >= 0.5","LKB > 1 Mt",
#                        "0.01< LHR <0.10","LHR >= 0.10","worst case"), cex=0.75)
axis(1,at=1:7,labels=c("1 (best)","7","13","2","3","5","18 (worst)"))
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)
## ADD BLUE LINE
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added

dev.copy(
  png,
  filename = file.path(out.dir, file = "./Watters_Penguin_performance.png"),
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

axis(1,at=1:18,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"),cex=0.75)
# axis(1, at=1:18, labels=FALSE)
# text(x=1:18, y=par("usr")[3]-0.3, labels=x_labels, srt=90, adj=1, xpd=TRUE, cex=0.65)

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
  filename = file.path(out.dir, file = "Watters_Penguin_performance_imputed.png"),
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
xx$impute.me<-ifelse(is.na(xx$survey) & xx$season=="S",1,0)
xx$imputed.hr<-hr.imputed.summ$statistics[,1]
xx$survey[xx$impute.me==1]<-xx$catch[xx$impute.me==1]/xx$imputed.hr[xx$impute.me==1]
xx$hr.class<-ifelse(xx$catch/xx$survey<=0.01,1,ifelse(xx$catch/xx$survey>=0.1,3,2))
xx$kb.class<-ifelse(xx$survey<=1000000,1,2)
xx$oni.class<-as.numeric(factor(xx$oni.class))
xx$case<-ifelse(xx$oni.class==1 & xx$kb.class==1 & xx$hr.class==1,1,
               ifelse(xx$oni.class==1 & xx$kb.class==2 & xx$hr.class==1,2,
               ifelse(xx$oni.class==1 & xx$kb.class==1 & xx$hr.class==2,3,
               ifelse(xx$oni.class==1 & xx$kb.class==2 & xx$hr.class==2,4,
               ifelse(xx$oni.class==1 & xx$kb.class==1 & xx$hr.class==3,5,
               ifelse(xx$oni.class==1 & xx$kb.class==2 & xx$hr.class==3,6,
               ifelse(xx$oni.class==2 & xx$kb.class==1 & xx$hr.class==1,7,
               ifelse(xx$oni.class==2 & xx$kb.class==2 & xx$hr.class==1,8,
               ifelse(xx$oni.class==2 & xx$kb.class==1 & xx$hr.class==2,9,
               ifelse(xx$oni.class==2 & xx$kb.class==2 & xx$hr.class==2,10,
               ifelse(xx$oni.class==2 & xx$kb.class==1 & xx$hr.class==3,11,
               ifelse(xx$oni.class==2 & xx$kb.class==2 & xx$hr.class==3,12,
               ifelse(xx$oni.class==3 & xx$kb.class==1 & xx$hr.class==1,13,
               ifelse(xx$oni.class==3 & xx$kb.class==2 & xx$hr.class==1,14,
               ifelse(xx$oni.class==3 & xx$kb.class==1 & xx$hr.class==2,15,
               ifelse(xx$oni.class==3 & xx$kb.class==2 & xx$hr.class==2,16,
               ifelse(xx$oni.class==3 & xx$kb.class==1 & xx$hr.class==3,17,18)))))))))))))))))


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
points(jitter(xx$case[xx$impute.me==1],amount=0.2),xx$index[xx$impute.me==1],cex=0.8,col="red")

abline(h=0)
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
## ADD BLUE LINE
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added

# Add custom x-axis with rotated labels
axis(1, at=1:18)
#text(x=1:18, y=par("usr")[3]-0.3, labels=x_labels, srt=90, adj=1, xpd=TRUE, cex=0.65)

# Add legend for ONI classes
legend("topright", 
       legend=c("ONI < -0.5C", "-0.5C < ONI < 0.5C", "ONI > 0.5C"),
       fill=c("lightblue", "white", "lightpink"),
       border="black",
       title="ONI Classes",
       cex=0.65)

dev.copy(
  png,
  filename = file.path(out.dir, file = "Watters_Posterior_predictive_dist_imputed.png"),
  width    = 8,
  height   = 6,
  units    = "in",
  res      = 300
)
dev.off()


##################################################################################################

###Recreation Summary Table 1, Watters et al. ---
#table original
d.org = read_csv("./A_Results/Watters_model/Imputed/hr_derived_s.csv",col_types="dddfd")

keepers = c(37:74)
pt = d.org %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarize(mean(value)) %>%
  dplyr::rename(.,c(mean = `mean(value)`)) %>%
  dplyr::slice(37:74)

best.case.names = c("Best case", "ONI Neutral", "ONI warm", "LKB High", "LHR medium", "LHR high", "Worst case")
best.case = c(NA,pt$mean[4], pt$mean[5], pt$mean[1], pt$mean[2], pt$mean[3], pt$mean[6])
oni.range = c(NA,NA,NA, pt$mean[9], pt$mean[7], pt$mean[8], NA)
oni.warm = c(NA, NA, NA, pt$mean[12], pt$mean[10], pt$mean[11], NA)
long.term.posterior = c(pt$mean[13], pt$mean[17], pt$mean[18], pt$mean[14], pt$mean[15], pt$mean[16], pt$mean[19])
long.term.predicted = c(pt$mean[32], pt$mean[36], pt$mean[37], pt$mean[33], pt$mean[34], pt$mean[35], pt$mean[38])

table.original.scenario = round(tibble(best.case, oni.range, 
                                       oni.warm, long.term.posterior,long.term.predicted),2) %>%
  dplyr::mutate(Effects = best.case.names) %>%
  dplyr::select(6,1,2,3,4,5)



##################################################################################################

# misc stuff

junk2<-make.localhr.data(plot.winter=TRUE)
junk2$impute.me<-ifelse(is.na(junk2$survey)&junk2$season=="S",1,0)
junk2$imputed<-exp(ifelse(junk2$impute.me==0,NA,ifelse(junk2$gSSMU==1&junk2$sam.sign=="Neg",hr.params.summ$statistics[7,1],
                                                   ifelse(junk2$gSSMU==2&junk2$sam.sign=="Neg",hr.params.summ$statistics[8,1],
                                                          ifelse(junk2$gSSMU==1&junk2$sam.sign=="Pos",hr.params.summ$statistics[9,1],
                                                                 hr.params.summ$statistics[10,1])))))

library(lattice)
# plot the time series
xyplot(index~cal.yr|season,data=junk2,horizontal=FALSE,aspect=0.25,panel=function(x,y,subscripts,...,Z=junk2$survey,IMP=junk2$imputed){
  z<-(Z[subscripts]-mean(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE))/sd(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE)
  imp<-(IMP[subscripts]-mean(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE))/sd(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE)
  panel.xyplot(x,y,...,pch=16,col="black")
  panel.points(as.numeric(x),z,col="red",pch=16,cex=1.25)
  panel.points(as.numeric(x),imp,col="red",pch=1,cex=1.25)
  panel.abline(h=0,lty=2)},
  ylim=c(-3.5,4.5),layout=c(1,2),xlab="year",ylab="std monitoring index")

# plot the time series with env indices (panel for Fig 1)
xyplot(index~cal.yr|season+PROJECT,data=junk2,horizontal=FALSE,aspect=0.5,panel=function(x,y,subscripts,...,Z=junk2$survey,IMP=junk2$imputed,ONI=junk2$oni,SAM=junk2$sam){
  z<-(Z[subscripts]-mean(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE))/sd(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE)
  imp<-(IMP[subscripts]-mean(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE))/sd(c(Z[subscripts],IMP[subscripts]),na.rm=TRUE)
  panel.xyplot(x,y,...,pch=16,col="black")
  panel.points(as.numeric(x),z,col="red",pch=16,cex=1.25)
  panel.points(as.numeric(x),imp,col="red",pch=1,cex=1.25)
  tt.oni<-data.frame(as.numeric(x),ONI[subscripts])
  tt.oni<-tt.oni[order(tt.oni[,1]),]
  panel.lines(tt.oni[,1],tt.oni[,2],col="blue")
  tt.sam<-data.frame(as.numeric(x),SAM[subscripts])
  tt.sam<-tt.sam[order(tt.sam[,1]),]
  panel.lines(tt.sam[,1],tt.sam[,2],col="dark green")
  panel.abline(h=0,lty=2)},
  ylim=c(-3.5,4.5),layout=c(2,2),xlab="year",ylab="index")


# some stuff to tabulate local harvest rates by year (including imputed estimates)

junk$impute.me<-ifelse(is.na(junk$survey)&junk$season=="S",1,0)
junk$survey<-ifelse(junk$impute.me==0,junk$survey,exp(ifelse(junk$gSSMU==1&junk$sam.sign=="Neg",hr.params.summ$statistics[7,1],
                                                       ifelse(junk$gSSMU==2&junk$sam.sign=="Neg",hr.params.summ$statistics[8,1],
                                                              ifelse(junk$gSSMU==1&junk$sam.sign=="Pos",hr.params.summ$statistics[9,1],
                                                                     hr.params.summ$statistics[10,1])))))
junk$hr<-junk$catch/junk$survey
junk$hihr<-(junk$hr>=0.10)

# table for supplementary info
jj<-unique(junk[,c(7,8,6,11,17,16)])
names(jj)<-c("calendar.year","stratum","season","catch","LHR","imputed")
jj$stratum<-ifelse(jj$stratum==1,"Bransfield","Drake")
jj$season<-ifelse(jj$season=="S","Summer","Winter")
jj$imputed<-ifelse(jj$imputed==1,"Yes","No")
# no catch data for 2016
jj<-jj[jj$calendar.year!=2016,]
write.csv(jj[order(jj[,1],jj[,2],jj[,3]),],file="hr.csv",row.names = FALSE)



# variation in catch by season and decade (panel for Fig 3)
tt<-read.csv("c1.csv",header=TRUE,stringsAsFactors = FALSE)
tt<-tt[is.element(tt$AssignedSSMU,c("APBSE","APBSW","APDPE","APDPW","APE","APEI","APPA","APW")),]
tt$FishingSeason<-ifelse(tt$Month==12,tt$CalendarYear+1,tt$CalendarYear)
tt$season<-ifelse(is.element(tt$Month,c(10:12,1:3)),"S","W")
tt$decade<-ifelse(tt$FishingSeason<1990,"before 1990",
                  ifelse(tt$FishingSeason>1989&tt$FishingSeason<2000,"1990-1999",
                         ifelse(tt$FishingSeason>1999&tt$FishingSeason<2010,"2000-2009","after 2009")))
tt$decade<-ordered(tt$decade,levels=c("before 1990","1990-1999","2000-2009","after 2009"))
tt<-tapply(tt$TotalCatch,list(tt$decade,tt$season),sum)
tt<-data.frame(catch=as.numeric(tt),decade=rep(dimnames(tt)[[1]],dim(tt)[2]),season=rep(dimnames(tt)[[2]],each=dim(tt)[1]))
tt$decade<-ordered(tt$decade,levels=c("before 1990","1990-1999","2000-2009","after 2009"))
#
barchart(I(catch/1000)~season|decade,data=tt,layout=c(4,1),aspect=1,xlab="Season",ylab="Total catch (1000 t)")
