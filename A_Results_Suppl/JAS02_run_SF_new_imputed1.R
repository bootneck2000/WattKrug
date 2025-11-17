#### MODEL EXECUTION
## WATTERS ET AL (2020)

## Conditions: 
# A. Survey data filtered by transect (nm) coverage 
# B. Missing summer LKB imputed, using new median estimates by gSSMU

library(tidyverse)
library(lattice)

# Jump to Step 3 if already ran steps 1 & 2

#### First Step, load variables & Update 'survey' data ####
setwd("C:/Users/javie/OneDrive/R-Git projects/WattKrug/")

# Same as JAS_run_Survey_filtered.R
junk <- readRDS("./A_Results_Suppl/Mod Survey filtered/junk.rds")

#### Second - Run Model ####
setwd("C:/Users/javie/OneDrive/R-Git projects/WattKrug/")
out.dir <- "./A_Results_Suppl/Mod2 LKB imputed1/"

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

## Here, modify hr.data to add new imputed data

hr.data<-list(index=as.vector(junk$index),
              survey=ifelse(is.na(junk$survey),1E12,junk$survey),
              catch=junk$catch,
              gssmu=junk$gSSMU,
              oniclass=as.numeric(factor(junk$oni.class)),
              samclass=as.numeric(factor(junk$sam.sign)),
              summer=junk$survey[junk$season=="S"],
              impute.me=ifelse(is.na(junk$survey)&junk$season=="S",1,0),
              meanlogsummer=tapply(log(junk$survey[junk$season=="S"]),
                                   list(junk$gSSMU[junk$season=="S"],junk$sam.sign[junk$season=="S"]),
                                   mean,na.rm=TRUE),
              sdlogsummer=sd(log(junk$survey[junk$season=="S"]),na.rm=TRUE),
              nobs=dim(junk)[1],
              nsummerobs=as.vector(table(junk$season)[1]),
              predX=pred.matrix)
## End

hr.params<-c("beta","mulogsummer","sigmalogsummer","sd.index","t.sd.index")

beta.init1<-rep(-1,6)
beta.init2<-rep(0,6)
beta.init3<-rep(1,6)


hr.inits<-list(list(beta=beta.init1,t.sd.index=0.1,.RNG.seed=123,
                    .RNG.name="base::Super-Duper"),
               list(beta=beta.init2,t.sd.index=1.0,.RNG.seed=456,
                    .RNG.name="base::Super-Duper"),
               list(beta=beta.init3,t.sd.index=1.9,.RNG.seed=789,
                    .RNG.name="base::Super-Duper"))

hr.derived<-c("index.new","mu.new","prob","prob.new")

hr.imputed<-"hr"

# now do the analysis
Sys.unsetenv("JAGS_HOME")  # let rjags auto-detect
library(coda)
library(rjags)

hr.jags<-jags.model(textConnection(modelstring),hr.data,hr.inits,n.chains=3,n.adapt=250)
# burn in for 150000 iterations
update(hr.jags, n.iter=50000)
hr.params.post<-coda.samples(hr.jags,hr.params,n.iter=10000,thin=25)
hr.derived.post<-coda.samples(hr.jags,hr.derived,n.iter=10000,thin=25)
hr.imputed.post<-coda.samples(hr.jags,hr.imputed,n.iter=10000,thin=25)
hr.params.summ<-summary(hr.params.post)
hr.derived.summ<-summary(hr.derived.post)
hr.imputed.summ<-summary(hr.imputed.post)

# write input/output
saveRDS(hr.params.post, file.path(out.dir, file = "hr.params.post.rds"))
saveRDS(hr.derived.post, file.path(out.dir, file = "hr.derived.post.rds"))
saveRDS(hr.imputed.post, file.path(out.dir, file = "hr.imputed.post.rds"))

saveRDS(hr.params.summ, file.path(out.dir, file = "hr.params.summ.rds"))
saveRDS(hr.derived.summ, file.path(out.dir, file = "hr.derived.summ.rds"))
saveRDS(hr.imputed.summ, file.path(out.dir, file = "hr.imputed.summ.rds"))

saveRDS(junk, file.path(out.dir, file = "junk.rds"))

# cat(capture.output(print(hr.params.post), file.path(out.dir, file="hr_params_post.txt")))
sink(file.path(out.dir, file="hr_params.txt"))
print(hr.params.post)
sink()
sink(file.path(out.dir, file="hr_derived.txt"))
print(hr.derived.summ)
sink()
sink(file.path(out.dir, file="hr_imputed.txt"))
print(hr.imputed.summ)
sink()


#### Step 3: PLOTTING  ####

out.dir <- "./A_Results_Suppl/Mod Survey filtered/"
setwd(out.dir)

# Read data files
hr.params.post <- readRDS("hr.params.post.rds")
hr.derived.post <- readRDS("hr.derived.post.rds")
hr.imputed.post2 <- readRDS("hr.imputed.post.rds")
junk <- readRDS("junk.rds")

hr.params.summ <- readRDS("hr.params.summ.rds")
hr.derived.summ <- readRDS("hr.derived.summ.rds")
hr.imputed.summ <- readRDS("hr.imputed.summ.rds")

# plot posterior expectations of marginal effects
require(ggmcmc)
hr.params.s<-ggs(hr.params.post)
hr.derived.s<-ggs(hr.derived.post)

# just want to copy hr.params.s to work with it for plotting diagnostics without screwing up the original object
# also get rid of chains for t.sd.index since this is not really a parameter of interest
HR.labels<-data.frame(Parameter=dimnames(hr.params.post[[1]])[[2]],
                      Label=c("alpha","beta[3]","beta[4]","beta[5]","beta[1]",
                              "beta[2]","K[B,-]","K[D,-]","K[B,+]","K[D,+]",
                              "sigma","phi","exclude"))

hr.params2.s<-ggs(hr.params.post,par_labels = HR.labels)
hr.params2.s<-hr.params2.s[hr.params2.s$ParameterOriginal!="t.sd.index",]

ggmcmc(hr.params.s, "diagnostics_hr_params_final.pdf")
ggmcmc(hr.derived.s, "diagnostics_hr_derived_final.pdf")

## Figure 2 paper Watters et al. (2020).

# reference (best case)
png("EP_summ_boxplot.png", width = 10, height = 7, units = "in", res = 300, pointsize = 10) # Added
par(mar = c(6, 5, 2, 1))  # optional margins # Added
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,7.5),ylab="Expected performance",
        xlab="",whisklty=1,boxwex=1,at=1)
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
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added blue line
dev.off() # Added


### Supplementary Figures

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

### Here, impute new values
xx$imputed.hr<-hr.imputed.summ$statistics[,1]
xx$survey[xx$impute.me==1]<-xx$catch[xx$impute.me==1]/xx$imputed.hr[xx$impute.me==1]

### End

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

#### Save xx
saveRDS(xx, "Data_Fig_S7.rds")  # save the data
####

# plot the data
png("std_performance_boxplot.png", width = 10, height = 7, units = "in", res = 300, pointsize = 10) # Added
par(mar = c(6, 5, 2, 1))  # optional margins # Added
plot(jitter(xx$case[xx$impute.me==0],amount=0.25),xx$index[xx$impute.me==0],type="n",xlim=c(0.65,18.35),
     ylim=c(-3.5,3.5),
     xlab="case",xaxt="n",ylab="std performance index",pch=16)
# add posterior predictive distributions as box plots
for(i in 1:18){
  boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(as.numeric(hr.derived.s$Parameter)==i),
          range=0,outline=FALSE,xaxt="n",whisklty=1,boxwex=1,at=i,add=TRUE)
}
# plot the data
points(jitter(xx$case[xx$impute.me==0&xx$season=="S"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="S"],
       cex=0.8,col="red",pch=16)
points(jitter(xx$case[xx$impute.me==0&xx$season=="W"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="W"],
       cex=0.8,col="blue",pch=16)
points(jitter(xx$case[xx$impute.me==1],amount=0.2),xx$index[xx$impute.me==1],cex=1,col="red")
axis(1,at=1:18)
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added
dev.off() # Added



### #  Added by Andy - plot all cases*

png("Expected_Performance_18.png", width = 10, height = 7, units = "in", res = 300, pointsize = 10) # Added
par(mar = c(6, 5, 2, 1))  # optional margins # Added

# # reference (best case)
boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,
        subset=(as.numeric(hr.derived.s$Parameter)==19),
        range=0,ylim=c(-2,2),xaxt="n",xlim=c(0.5,18.5),
        ylab="Expected performance",xlab="Case scenario",
        whisklty=1,boxwex=1,at=1)

boxplot(value~I(as.numeric(Parameter)),data=hr.derived.s,subset=(is.element(as.numeric(hr.derived.s$Parameter),c(20:36))),
        range=0,xaxt="n",yaxt="n",whisklty=1,boxwex=0.5,add=TRUE,at=2:18,col="gray80")

# points added by Javier (from  S7)
points(jitter(xx$case[xx$impute.me==0&xx$season=="S"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="S"],
       cex=0.8,col="red",pch=16)
points(jitter(xx$case[xx$impute.me==0&xx$season=="W"],amount=0.2),xx$index[xx$impute.me==0&xx$season=="W"],
       cex=0.8,col="blue",pch=16)
points(jitter(xx$case[xx$impute.me==1],amount=0.2),xx$index[xx$impute.me==1],cex=1,col="red")
##
axis(1,at=1:18,labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"))
abline(h=mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),lty=2)
abline(h=0)
abline(h=-mean(hr.derived.s$value[as.numeric(hr.derived.s$Parameter)==19]),col="blue",lty=2)  # Added
dev.off() # Added





#### misc stuff ####

junk2<-make.localhr.data(plot.winter=TRUE)
junk2$impute.me<-ifelse(is.na(junk2$survey)&junk2$season=="S",1,0)
junk2$imputed<-exp(ifelse(junk2$impute.me==0,NA,ifelse(junk2$gSSMU==1&junk2$sam.sign=="Neg",
                                                       hr.params.summ$statistics[7,1],
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
xyplot(index~cal.yr|season+PROJECT,data=junk2,horizontal=FALSE,aspect=0.5,panel=function(x,y,subscripts,...,Z=junk2$survey,
                                                                                         IMP=junk2$imputed,ONI=junk2$oni,SAM=junk2$sam){
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
junk$survey<-ifelse(junk$impute.me==0,junk$survey,exp(ifelse(junk$gSSMU==1&junk$sam.sign=="Neg",
                                                             hr.params.summ$statistics[7,1],
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
tt<-data.frame(catch=as.numeric(tt),decade=rep(dimnames(tt)[[1]],dim(tt)[2]),season=rep(dimnames(tt)[[2]],
                                                                                        each=dim(tt)[1]))
tt$decade<-ordered(tt$decade,levels=c("before 1990","1990-1999","2000-2009","after 2009"))
#
barchart(I(catch/1000)~season|decade,data=tt,layout=c(4,1),aspect=1,xlab="Season",ylab="Total catch (1000 t)")

           

