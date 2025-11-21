### Watter et al. 2020

# Model imputed LKB
# Re-creation of equations (1) and (2): mean(k) and LKBij

library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(scales)

#### First, import data from the model ####

# krill survey biomass
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)

# Filter low transect coverage (<10th percentil)
survey_filtered <- survey %>%
  filter(gSSMU %in% c(1, 2))
percentiles <- survey_filtered %>%
  group_by(gSSMU) %>%
  summarise(p10 = quantile(nmi.count, probs = 0.10, na.rm = TRUE))
percentiles$p10 <- ceiling(percentiles$p10/10)*10
survey_low_nmi <- survey_filtered %>%
  left_join(percentiles, by = "gSSMU") %>%
  filter(nmi.count >= p10)
survey <- survey_low_nmi
rm(percentiles, survey_filtered, survey_low_nmi)
### END changes

survey<-tapply(survey$biomass,list(survey$Year,survey$gSSMU),mean,na.rm=TRUE)
survey<-data.frame(cal.yr=rep(dimnames(survey)[[1]],dim(survey)[2]),
                   gSSMU=rep(dimnames(survey)[[2]],each=dim(survey)[1]),
                   survey=c(survey),stringsAsFactors = FALSE)
survey$season<-ifelse(survey$cal.yr<2012,"S","W")
survey$matchme<-paste(survey$cal.yr,survey$season,survey$gSSMU,sep="|")

## Adding Chinese survey data from Wang et al. 2025. ICES J Mar Sci, 82(8)
# Create dataframe 'Chinese_north' (gSSMU2)
Chinese_SSIW <- data.frame(
  cal.yr  = c("2013", "2015", "2016"),
  gSSMU   = rep("2", 3),
  survey    = c(418000, 400000, 896000),
  season  = rep("S", 3),
  stringsAsFactors = FALSE)

Chinese_BS <- data.frame(
  cal.yr  = c("2013", "2015", "2016"),
  gSSMU   = rep("1", 3),
  survey    = c(563000, 395000, 1035000),
  season  = rep("S", 3),
  stringsAsFactors = FALSE)

Chinese_survey <- rbind(Chinese_SSIW, Chinese_BS) 
Chinese_survey$matchme <- paste(Chinese_survey$cal.yr, Chinese_survey$season,
                                Chinese_survey$gSSMU,sep="|")
### Add Chinese survey data into 'survey'
survey <- rbind(survey, Chinese_survey) 

rm(Chinese_BS, Chinese_SSIW, Chinese_survey)
### END changes

survey <- survey %>% mutate(cal.yr = as.integer(cal.yr))
# survey <- survey %>% filter(season == "S")   # only imputed summer data
survey <- survey %>% mutate(gSSMU = as.integer(gSSMU))


# Import SAM data
sam<-read.csv("./Supplementary Files/sam.csv")
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$cal.yr<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$cal.yr,sam$season),mean)
sam<-data.frame(cal.yr=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))

sam <- sam %>% mutate(cal.yr = as.integer(cal.yr))
sam <- sam %>% filter(season == "S")   # only imputed summer data


# krill fishery catches
fishery<-read.csv("./Supplementary Files/c1.csv",header=TRUE,stringsAsFactors = FALSE)
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

fishery <- fishery %>% mutate(cal.yr = as.integer(cal.yr)) 
fishery <- fishery %>% filter(season == "S")   # only imputed summer data

# junk <- readRDS("./javier_analysis/Watters_model_output1/junk.rds")
# data <- junk %>% filter(season == "S") %>% filter(!is.na(survey)) %>% distinct(survey)

#### Imputing LKB data based on SAM ####

# join variables
data <- left_join(fishery, survey, by = c("cal.yr", "gSSMU", "season", "matchme")) 
data <- data %>% filter(gSSMU %in% c(1,2))
data <- left_join(data, sam, by = c("cal.yr", "season")) 
data <- data %>% mutate(sam.sign = (ifelse(sam<0, "Neg", "Pos")))

# Estimating mean(logLKB) for survey data
data$LnSurvey <- log(data$survey)

meanLogS <- data.frame(SAM.neg = numeric(2), SAM.pos = numeric(2))
sdLogS <- data.frame(SAM.neg = numeric(2), SAM.pos = numeric(2))

data.pos <- data %>% filter(sam.sign == "Pos")
data.neg <- data %>% filter(sam.sign == "Neg")

for(i in 1:2) {
  gS <- data.neg %>% filter(!is.na(survey)) %>% 
    filter(gSSMU == i) %>% filter(season == "S") %>%
    filter(cal.yr<2012)   # Estimate Mean/SD only with original data
  meanLogS[i,"SAM.neg"] <- mean(gS$LnSurvey)
  sdLogS[i, "SAM.neg"] <- sd(gS$LnSurvey)
}
for(i in 1:2) {
  gS <- data.pos %>% filter(!is.na(survey)) %>% 
    filter(gSSMU == i) %>% filter(season == "S") %>%
    filter(cal.yr<2012)   # Estimate Mean/SD only with original data
  meanLogS[i, "SAM.pos"] <- mean(gS$LnSurvey)
  sdLogS[i,"SAM.pos"] <- sd(gS$LnSurvey)
}
rm(gs)

# Imputing missing data
data <- data %>%
  mutate(
    # 0 if survey available, 1 if missing
    # impute.me = 1 # if_else(is.na(survey), 1, 0),  #Impute the whole series, for comparison
    
    # fill log survey
    LogSurvey_imputed = case_when(
      gSSMU == 1  & sam.sign == "Neg" ~ meanLogS[1, "SAM.neg"],  # remove from begining: 'impute.me == 1 &' 
      gSSMU == 1  & sam.sign == "Pos" ~ meanLogS[1, "SAM.pos"],
      gSSMU == 2  & sam.sign == "Neg" ~ meanLogS[2, "SAM.neg"],
      gSSMU == 2  & sam.sign == "Pos" ~ meanLogS[2, "SAM.pos"],
      TRUE ~ NA_real_
    ),
    
    # back-transform
    survey_imputed = exp(LogSurvey_imputed)
  )

#### Statistics and Plots  ####

### Statistical Analysis 
library(coin)

g1 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 1) %>% filter(cal.yr<2012)
g1 <- g1 %>% mutate(sam.sign = as.factor(sam.sign))
g1.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g1, distribution = "approximate")

g2 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 2) %>% filter(cal.yr<2012)
g2 <- g2 %>% mutate(sam.sign = as.factor(sam.sign))
g2.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g2, distribution = "approximate")

#### Table 1 

Table.01 <- tibble(
  "Area" = c("gSSMU1", "gSSMU1", "gSSMU2", "gSSMU2"),
  "SAM" = c("Negative", "Positive", "Negative", "Positive"),
  "Mean(LKB)" = NA_real_,
  "SD(LKB)"   = NA_real_,
  "p-value" = NA_real_
)

Table.01$`Mean(LKB)`[1] <- exp(meanLogS[1,1])
Table.01$`Mean(LKB)`[2] <- exp(meanLogS[1,2])
Table.01$`Mean(LKB)`[3] <- exp(meanLogS[2,1])
Table.01$`Mean(LKB)`[4] <- exp(meanLogS[2,2])

Table.01$`SD(LKB)`[1] <- exp(sdLogS[1,1])
Table.01$`SD(LKB)`[2] <- exp(sdLogS[1,2])
Table.01$`SD(LKB)`[3] <- exp(sdLogS[2,1])
Table.01$`SD(LKB)`[4] <- exp(sdLogS[2,2])

Table.01$`p-value`[1] <- pvalue(g1.stat)
Table.01$`p-value`[3] <- pvalue(g2.stat)

knitr::kable(Table.01, digits = 3)


### Plots



### Fig S2: survey LKB vs. imputed data gSSMU1 and 2
Chinese <- survey %>% filter(cal.yr > 2012) %>% filter(season == "S") 
FigS2_1 <- survey %>%
  filter(gSSMU == "1") %>%
  ggplot(aes(x = cal.yr, y = survey, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(data = filter(survey, gSSMU == "1", cal.yr <= 2012),
             shape = 16, size = 5) +   # filled circles
  geom_point(data = filter(survey, gSSMU == "1", season == "W"),
             shape = 16, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) +
  geom_point(data = filter(data, gSSMU == "1"),
             aes(x = cal.yr, y = survey_imputed, color = season, group = season),
             shape = 21, size = 3.5, na.rm = TRUE) +
  geom_point(data = filter(Chinese, gSSMU == "1"),
             shape = 17, size = 5, color = "red") +   # filled triangles
  scale_y_continuous(labels = scales::label_comma(),
                     breaks = seq(0, max(survey$survey, na.rm = TRUE), by = 1e6)) +
  labs(title = "Krill Survey - BS",
       x = "Year", y = "Biomass (tons)", color = "Season") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); FigS2_1


FigS2_2 <- survey %>%
  filter(gSSMU == "2") %>%
  ggplot(aes(x = cal.yr, y = survey, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(data = filter(survey, gSSMU == "2", cal.yr <= 2012),
             shape = 16, size = 5) +   # filled circles
  geom_point(data = filter(survey, gSSMU == "2", season == "W"),
             shape = 16, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) +
  geom_point(data = filter(data, gSSMU == "2"),
             aes(x = cal.yr, y = survey_imputed, color = season, group = season),
             shape = 21, size = 3.5, na.rm = TRUE) +
  geom_point(data = filter(Chinese, gSSMU == "2"),
             shape = 17, size = 5, color = "red") +   # filled triangles
  scale_y_continuous(labels = scales::label_comma(),
                     breaks = seq(0, max(survey$survey, na.rm = TRUE), by = 1e6)) +
  labs(title = "Krill Survey - SSWI-EI",
       x = "Year", y = "Biomass (tons)", color = "Season") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); FigS2_2



### Estimating Frequency distributions for LKB and imputed ####

# Survey data
LKB_survey <- data 
params_ln <- data %>%
  filter(!is.na(LnSurvey),
         sam.sign %in% c("Neg", "Pos"),
         gSSMU %in% c(1, 2)) %>%
  group_by(sam.sign, gSSMU) %>%
  summarise(
    n       = n(),
    mu_hat  = mean(LnSurvey),
    sd_hat  = sd(LnSurvey),
    .groups = "drop"
  )
params_ln

# 1. Data for histograms
hist_data <- data %>%
  filter(!is.na(LnSurvey),
         sam.sign %in% c("Neg", "Pos"),
         gSSMU %in% c(1, 2)) %>%
  mutate(
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU)
  )

# 2. Data for fitted normal curves, built only from params_ln
curve_data <- params_ln %>%
  mutate(
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU)
  ) %>%
  rowwise() %>%
  mutate(
    x = list(seq(mu_hat - 4 * sd_hat, mu_hat + 4 * sd_hat, length.out = 200))
  ) %>%
  unnest(x) %>%
  mutate(
    density = dnorm(x, mean = mu_hat, sd = sd_hat)
  ) %>%
  ungroup()

# 3. Plot: histogram of data + normal line from params_ln
x_min <- min(hist_data$LnSurvey, na.rm = TRUE)
x_max <- max(hist_data$LnSurvey, na.rm = TRUE)
survey_hist <- ggplot(hist_data, aes(x = LnSurvey)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.5,
                 color = "black",
                 fill  = "grey80") +
  geom_line(data = curve_data,
            aes(x = x, y = density),
            linewidth = 1) +
  facet_grid(sam.sign ~ gSSMU) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  labs(x = "LnSurvey", y = "Density") +
  theme_minimal()
survey_hist



 
## Figure S3: 
ggplot(data) +
  geom_density(aes(x = log(survey), color = "Original", group = sam.sign), linewidth = 1.2, alpha = 0.6) +
  geom_density(aes(x = log(survey_imputed), color = "Imputed", group = sam.sign), linewidth = 1.2, alpha = 0.6) +
  facet_wrap(~ gSSMU, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Original" = "steelblue", "Imputed" = "orange")) +
  labs(title = "Survey vs Imputed survey distributions by gSSMU",
       x = "Survey value",
       y = "Density",
       color = "Dataset") +
  theme_minimal(base_size = 14)



## ---- Priors -----------------------------------------------------------
# # mulogsummer[i,j] ~ dunif(0.1 * meanlogsummer[i,j], 10 * meanlogsummer[i,j])
# mulogsummer <- matrix(NA_real_, nrow = 2, ncol = 2)
# for (i in 1:2) {
#   for (j in 1:2) {
#     mulogsummer[i, j] <- runif(
#       n   = 1,
#       min = 0.1 * meanLogS[i, j],
#       max = 10  * meanLogS[i, j]
#     )
#   }
# }
# 
# # sigmalogsummer ~ dunif(0.1 * sdlogsummer, 10 * sdlogsummer)
# sigmalogsummer <- runif(
#   n   = 1,
#   min = 0.1 * sdLogS,
#   max = 10  * sdLogS
# )
# taulogsummer <- sigmalogsummer^(-2)
# 
# 
# 
# for(i in nrow(missLKB)) {
#   
# }
# 
# 