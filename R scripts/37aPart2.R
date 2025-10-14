junk<-make.localhr.data()
# write out data
write.csv(junk, file = "junk.csv")  
modelstring<-"
 
      model{
    for(i in 1:nsummerobs){
      lower[i]<-max(10000,catch[i])
      summer[i]~dlnorm(mulogsummer[gssmu[i],samclass[i]],taulogsummer) T(lower[i],100000000)
    }
    
    for(i in 1:2){ # two gSSMUs #change here
      for(j in 1:2){ # two SAM classes
        mulogsummer[i,j]~dunif(0.1*meanlogsummer[i,j],10*meanlogsummer[i,j])
      }
    }
    taulogsummer<-pow(sigmalogsummer,-2)
    sigmalogsummer~dunif(0.1*sdlogsummer,10*sdlogsummer)
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

   
    for(i in 1:nobs){
      X[i,1]<-1.0     # intercept
      X[i,2]<-equals(bclass[i],2)-equals(bclass[i],1) # b2
      X[i,3]<-equals(hrclass[i],2)-equals(hrclass[i],1) # hr2
      X[i,4]<-equals(hrclass[i],3)-equals(hrclass[i],1) # hr3
      X[i,5]<-equals(oniclass[i],2)-equals(oniclass[i],1) # o2
      X[i,6]<-equals(oniclass[i],3)-equals(oniclass[i],1) # o3
    }
  
 
    for(i in 1:nobs){
      index[i]~dnorm(mu[i],tau.index)
      mu[i] <- inprod(X[i,],beta[])
    }
    
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
