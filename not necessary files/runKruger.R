library(tidyverse)
library(rgdal)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read penguin data:
## Downloaded directly:
peng <- read.csv('./Kruger/selected_counts_V_3.0.csv', header=T)

## Subset according to kruger et al.:
peng <- peng %>% filter(common.name!='adelie penguin', 
                        count.type=='nests', 
                        (month=='11' | month=='12'),
                        as.numeric(year)>=1980, as.numeric(year)<=2017)
names(peng) <- gsub('.EPSG.4326', '', names(peng))

peng$date <- as.POSIXct(strptime(paste(peng$year, peng$month, peng$day), '%Y%m%d'), tz='UTC')
peng$count <- as.numeric(peng$count)

## From Kruger et al. supplementary:
krug <- tabulizer::extract_tables('./Kruger/13280_2020_1386_MOESM1_ESM.pdf', pages=c(4:13))
krug[[1]] <- krug[[1]][-c(1:2),]
krug[[1]] <- t(apply(krug[[1]], 1, function(x) unlist(strsplit(x, ' '))))
krug <- as.data.frame(do.call('rbind', krug))
names(krug) <- c('common.name', 'site.name', 'Longitude', 'Latitude', 'year', 'catch', 'SAM', 'count', 'lambda') 
for(i in 3:length(krug)) krug[,i] <- as.numeric(krug[,i])
krug$binLambda <- -sign(krug$lambda)
krug$binLambda[which(krug$binLambda<0)] <- 0
## library(lmerTest)

krug$binL <- krug$binLambda


library(lme4)
krugGLMERchin <- glmer(binLambda~log(catch+1)*SAM + (1|site.name), data=krug %>% filter(common.name=='Chinstrap'), family=binomial)
krugGLMERgent <- glmer(binLambda~log(catch+1)*SAM + (1|site.name), data=krug %>% filter(common.name=='Gentoo'), family=binomial)

library(lmerTest)
krugLMERchin <- lmer(binLambda~log(catch+1)*SAM + (1|site.name), data=krug %>% filter(common.name=='Chinstrap'))
krugLMERgent <- lmer(binLambda~log(catch+1)*SAM + (1|site.name), data=krug %>% filter(common.name=='Gentoo'))


sjPlot::plot_model(krugLMERchin, type='int')
sjPlot::plot_model(krugLMERgent, type='int')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create spatial object:
peng.sp <- peng
coordinates(peng.sp) <- c('Longitude','Latitude')
proj4string(peng.sp) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')

## Get CCAMLR Small-scale Management Units
## by overlaying peng.sp with ssmu:

ssmu <- readOGR('./ssmu-shapefile-WGS84/ssmu-shapefile-WGS84.shp')
peng$ssmu <- over(peng.sp, ssmu)$Name
peng$ssmu <- gsub("SSMU 48.1 ", "", peng$ssmu)

for(i in which(is.na(peng$ssmu))) {
  closest <- which.min(geosphere::distHaversine(peng.sp[i,], 
                                     peng.sp[which(!is.na(peng$ssmu)),]))
  peng$ssmu[i] <- peng$ssmu[which(!is.na(peng$ssmu))][closest]
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Summarize by site and year:

peng.s <- peng %>% 
  group_by(ssmu, common.name, site.name, year) %>% 
  summarize(mean=mean(count), max=max(count), min=min(count)) %>% 
  arrange(common.name, site.name, year)

## Calculate growth rate
peng.s$lambda <- rep(NA, nrow(peng.s))

for(i in 2:nrow(peng.s)) {
  peng.s$lambda[i] <- ((peng.s$mean[i]/peng.s$mean[i-1])/(peng.s$year[i]-peng.s$year[i-1]))-1
}

## Binary pop change!!

peng.s$binLambda <- c(NA, sign(diff(peng.s$mean)))

peng.s$lambda[match(unique(peng.s$site.name), peng.s$site.name)] <- NA
peng.s$binLambda[match(unique(peng.s$site.name), peng.s$site.name)] <- NA


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get environmental data:

sam <- readLines('http://www.nerc-bas.ac.uk/public/icd/gjma/newsam.1957.2007.txt')[-c(1,2)]
sam <- as.data.frame(do.call('rbind', lapply(sam, function(x) {
  spl <- unlist(strsplit(x, ' '))
  spl <- spl[which(nchar(spl)>0)]
  if(length(spl)<13) spl <- c(spl, rep(NA, 13-length(spl)))
  as.numeric(spl)
  })))
names(sam) <- c('year', )