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

# *******************************************************************************************************
# *************************************** Plot the input ************************************************
# *******************************************************************************************************

plot(as.vector(junk$index))
plot(ifelse(is.na(junk$survey),1E12,junk$survey))
plot(junk$catch)
plot(junk$gSSMU)
plot(as.numeric(factor(junk$oni.class)))
plot(as.numeric(factor(junk$sam.sign)))
plot(junk$survey[junk$season=="S"])
plot(ifelse(is.na(junk$survey)&junk$season=="S",1,0))
plot(tapply(log(junk$survey[junk$season=="S"]),list(junk$gSSMU[junk$season=="S"],junk$sam.sign[junk$season=="S"]),mean,na.rm=TRUE))
plot(sd(log(junk$survey[junk$season=="S"]),na.rm=TRUE))
plot(dim(junk)[1])
plot(as.vector(table(junk$season)[1]))

# *******************************************************************************************************

hr.params<-c("beta","mulogsummer","sigmalogsummer","sd.index","t.sd.index")

beta.init1<-rep(-1,6)
beta.init2<-rep(0,6)
beta.init3<-rep(1,6)


hr.inits<-list(list(beta=beta.init1,t.sd.index=0.1,.RNG.seed=123,.RNG.name="base::Super-Duper"),
               list(beta=beta.init2,t.sd.index=1.0,.RNG.seed=456,.RNG.name="base::Super-Duper"),
               list(beta=beta.init3,t.sd.index=1.9,.RNG.seed=789,.RNG.name="base::Super-Duper"))

hr.derived<-c("index.new","mu.new","prob","prob.new")

hr.imputed<-"hr"

# write out the input (to check how we are doing)
write.csv(hr.data$index, file = "index.csv")
write.csv(hr.data$survey, file = "survey.csv")
write.csv(hr.data$catch, file = "catch.csv")
write.csv(hr.data$gssmu, file = "gssmu.csv")
write.csv(hr.data$oniclass, file = "oniclass.csv")
write.csv(hr.data$samclass, file = "samclass.csv")
write.csv(hr.data$summer, file = "summer.csv")
write.csv(hr.data$impute.me, file = "imputeme.csv")

write.csv(hr.data$meanlogsummer, file = "meanlogsummer.csv")
write.csv(hr.data$sdlogsummer, file = "sdlogsummer.csv")
write.csv(hr.data$nobs, file = "nobs.csv")
write.csv(hr.data$nsummerobs, file = "nsummerobs.csv")
write.csv(hr.data$predX, file = "predX.csv")


# now do the analysis
library(coda)
library(rjags)

hr.jags<-jags.model(textConnection(modelstring),hr.data,hr.inits,n.chains=3,n.adapt=250000)
# burn in for 150000 iterations
update(hr.jags, n.iter=500000)
hr.params.post<-coda.samples(hr.jags,hr.params,n.iter=125000,thin=25)
hr.derived.post<-coda.samples(hr.jags,hr.derived,n.iter=125000,thin=25)
hr.imputed.post<-coda.samples(hr.jags,hr.imputed,n.iter=125000,thin=25)
hr.params.summ<-summary(hr.params.post)
hr.derived.summ<-summary(hr.derived.post)
hr.imputed.summ<-summary(hr.imputed.post)

# write input/output
# cat(capture.output(print(hr.params.post), file="hr_params_post.txt"))
sink("hr_params.txt")
print(hr.params.post)
sink()
sink("hr_derived.txt")
print(hr.derived.summ)
sink()
sink("hr_imputed.txt")
print(hr.imputed.summ)
sink()


require(ggmcmc)
hr.params.s<-ggs(hr.params.post)
hr.derived.s<-ggs(hr.derived.post)

# just want to copy hr.params.s to work with it for plotting diagnostics without screwing up the original object
# also get rid of chains for t.sd.index since this is not really a parameter of interest
HR.labels<-data.frame(Parameter=dimnames(hr.params.post[[1]])[[2]],
                      Label=c("alpha","beta[3]","beta[4]","beta[5]","beta[1]","beta[2]","K[B,-]","K[D,-]","K[B,+]","K[D,+]","sigma","phi","exclude"))

hr.params2.s<-ggs(hr.params.post,par_labels = HR.labels)
hr.params2.s<-hr.params2.s[hr.params2.s$ParameterOriginal!="t.sd.index",]

ggmcmc(hr.params.s,file="diagnostics_hr_params_final.pdf")
ggmcmc(hr.derived.s,file="diagnostics_hr_derived_final.pdf")


# ************************************************************************************

# plot posterior expectations of marginal effects

# ************************************************************************************

# Figure 2

# reference (best case)
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,7.5),ylab="expected performance",whisklty=1,boxwex=1,at=1)
# ONI
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(25,31))),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=2:3,col="gray80")
# biomass
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==20),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=1,add=TRUE,at=4,col="gray40",medcol="white")
# harvest rate
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(21,23))),
        range=0,yaxt="n",xaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=5:6,col="black",medcol="white")
# worst case with PARAMETER 36 chosen - so with "warm" ONI.  Should be (hr.derived.s$Parameter)==12
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==36),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=1,add=TRUE,at=7)
axis(1,at=1:7,labels=c("reference","-0.5 < ONI < 0.5","ONI >= 0.5","LKB > 1 Mt","0.01 < LHR < 0.10","LHR >= 0.10","worst case"))
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)

write.csv(hr.derived.s, file = "hr_derived_s.csv")  
# # ***********************************************************************************************
# # *************************************** Added by me - plot all cases***************************
# # ***********************************************************************************************
# # reference (best case)
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,18.5),ylab="expected performance",whisklty=1,boxwex=1,at=1)

boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(20:36))),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=2:18,col="gray80")

axis(1,at=1:18,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)

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
# plot the data
plot(jitter(xx$case[xx$impute.me==0],amount=0.25),xx$index[xx$impute.me==0],type="n",xlim=c(0.65,18.35),ylim=c(-3.5,3.5),xlab="case",xaxt="n",ylab="std performance index",pch=16)
# add posterior predictive distributions as box plots
for(i in 1:18){
  boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==i),
          range=1.5,outline=FALSE,xaxt="n",whisklty=1,boxwex=1,at=i,add=TRUE)
}
# plot the data
points(jitter(xx$case[xx$impute.me==0&xx$season=="S"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="S"],cex=0.5,col="red",pch=16)
points(jitter(xx$case[xx$impute.me==0&xx$season=="W"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="W"],cex=0.5,col="blue",pch=16)
points(jitter(xx$case[xx$impute.me==1],amount=0.2),xx$index[xx$impute.me==1],cex=0.5,col="red")
axis(1,at=1:18)
abline(h=0)






##################################################################################################
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


# some stuff to tablulate local harvest rates by year (including imputed estimates)

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
