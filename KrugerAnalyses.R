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

peng <- peng %>% filter(count.type=='nests', 
                        (month=='11' | month=='12'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read the data in the kruger et al. supplement:

krug <- tabulizer::extract_tables('./Kruger/13280_2020_1386_MOESM1_ESM.pdf', pages=c(4:13))
krug[[1]] <- krug[[1]][-c(1:2),]
krug[[1]] <- t(apply(krug[[1]], 1, function(x) unlist(strsplit(x, ' '))))
krug <- as.data.frame(do.call('rbind', krug))
names(krug) <- c('common.name', 'site.name', 'Longitude', 'Latitude', 'year', 'catch', 'SAM', 'count', 'lambda') 
for(i in 3:length(krug)) krug[,i] <- as.numeric(krug[,i])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare the downloaded and supplementary data:
compMat <- cbind(as.vector(by(peng, peng$common.name, function(y) length(unique(y$site.name)))),
                 as.vector(by(krug, krug$common.name, function(y) length(unique(y$site.name)))))
dimnames(compMat) <- list(Species=c('Chinstrap', 'Gentoo'),
                          Dataset=c('MAAPD download', 'Kruger supplementary'))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Check distribution of intervals between surveys:
krug <- krug[order(krug$site.name, krug$year),]
gaps <- by(krug, list(krug$site.name, krug$common.name), function(x) {
  yr.diff <- diff(x$year)
  n.per <- length(yr.diff)
  n.gap <- length(which(yr.diff>1))
  mean.diff <- mean(yr.diff[which(yr.diff>1)])
##  c(n.per, n.gap, mean.diff)
  yr.diff
})

gaps <- krug %>% 
  group_by(site.name, common.name) %>%
  summarize(yrdiff=diff(year)) %>% as.data.frame()

gaptab <- table(gaps$yrdiff, gaps$common.name)
gaptab <- round(100*gaptab/apply(gaptab, 2, sum), 2)

for(i in 2:dim(gaptab)[1]) {
  gaptab[i,1] <- sum(gaptab[c(i:dim(gaptab)[1]),1])
  gaptab[i,2] <- sum(gaptab[c(i:dim(gaptab)[1]),2])
}

rownames(gaptab)[-c(1,2)] <- paste('>', as.numeric(rownames(gaptab)[-c(1,2)])-1, ' years', sep='') 
rownames(gaptab)[1] <- '1 year'
rownames(gaptab)[2] <- '>1 year'

