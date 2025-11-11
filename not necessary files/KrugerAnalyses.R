library(tidyverse)
library(lme4)
library(lmerTest)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read MAPPPD data:
peng <- read.csv('./Kruger/AllData_V_3.0.csv', header=T)
names(peng) <- gsub('_epsg_4326', '', names(peng))
names(peng) <- gsub('_', '.', names(peng), fixed=T)

## Subset according to kruger et al.:
peng <- peng %>% filter((common.name=='chinstrap penguin' | common.name=='gentoo penguin'), as.numeric(year)>=1980, as.numeric(year)<=2017,
                        latitude>(-66), longitude>-66, longitude<(-54))

peng$date <- as.POSIXct(strptime(paste(peng$year, peng$month, peng$day), '%Y%m%d'), tz='UTC')
peng$count <- as.numeric(peng$count)

## Save initial dataset comparison table:

peng.init <- peng

peng.n <- peng %>% 
  group_by(site.name) %>% 
  count() %>% filter(n>1)

peng <- right_join(peng, peng.n, by='site.name') 

peng <- peng %>% filter(count.type!='chicks', 
                        (month=='11' | month=='12'))

names(peng)[grep('penguin.count', names(peng))] <- 'count'

peng$count[which(peng$count.type=='adults')]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read the data in the kruger et al. supplement:
load("./Kruger/kruger_tables_extracted.RData") #tabulapdf::extract_area followed by rbdind

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare the downloaded and supplementary data:
compMat <- cbind(as.vector(by(peng, peng$common.name, function(y) length(unique(y$site.name)))),
                 as.vector(by(krug, krug$species, function(y) length(unique(y$loc)))))
dimnames(compMat) <- list(Species=c('Chinstrap', 'Gentoo'),
                          Dataset=c('MAAPD download', 'Kruger supplementary'))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Check distribution of intervals between surveys:
krug <- krug[order(krug$loc, krug$year),]
gaps <- by(krug, list(krug$loc, krug$species), function(x) {
  yr.diff <- diff(x$year)
  n.per <- length(yr.diff)
  n.gap <- length(which(yr.diff>1))
  mean.diff <- mean(yr.diff[which(yr.diff>1)])
##  c(n.per, n.gap, mean.diff)
  yr.diff
})

gaps <- krug %>% 
  group_by(loc, species) %>%
  summarize(yrdiff=diff(year)) %>% as.data.frame()

gaptab <- table(gaps$yrdiff, gaps$species)
gaptab <- round(100*gaptab/apply(gaptab, 2, sum), 2)

for(i in 2:dim(gaptab)[1]) {
  gaptab[i,1] <- sum(gaptab[c(i:dim(gaptab)[1]),1])
  gaptab[i,2] <- sum(gaptab[c(i:dim(gaptab)[1]),2])
}

rownames(gaptab)[-c(1,2)] <- paste('>', as.numeric(rownames(gaptab)[-c(1,2)])-1, ' years', sep='') 
rownames(gaptab)[1] <- '1 year'
rownames(gaptab)[2] <- '>1 year'

max.krug<- krug %>% filter(loc==gaps$loc[which.max(gaps$yrdiff)], species==gaps$species[which.max(gaps$yrdiff)])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create binary population change index:
krug$binLambda <- -sign(krug$lambda)
krug$binLambda[which(krug$binLambda<0)] <- 0


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Check frequency of zero catches:
noCatch <- inner_join(krug %>% group_by(species) %>% summarise(len=length(unique(loc))),krug %>% filter(catch>0) %>% group_by(species) %>% summarise(len=length(unique(loc))), by='species')
names(noCatch) <- c('Species', 'All colonies', 'Colonies with non-zero catches')


