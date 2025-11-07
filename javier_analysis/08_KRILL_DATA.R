##### KRILL SURVEY AND FISHERY DATA

## Watters Code:

# krill survey biomass
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)
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
#print(str(fishery))
#


# ==================== ANALYSIS  ====================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})
suppressWarnings(RNGkind("default"))  # silence harmless RNG warnings on some R builds

# -------------------- Preconditions & setup ------------------------------
# stopifnot(exists("survey"), is.data.frame(survey))
stopifnot(exists("fishery"), is.data.frame(fishery))

# --- Output locations ---------------------------
out_dir <- "./javier_analysis/Krill_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------ Fishing ratio (Catch proportion of Krill Biomass) ---- 

# --- Clean & aggregate survey ---
source("./javier_analysis/survey_filtered.R") # estimates 'survey' after filtering poor coverage

# --- Clean & aggregate fishery ---
# Run after Watters code: # krill fishery catches
fishery_by <- fishery %>%
  mutate(
    cal.yr = as.integer(cal.yr),          # fishery cal.yr may be num -> ensure integer
    gSSMU  = as.character(gSSMU),         # fishery gSSMU may be num -> make chr to join
    season = toupper(season)
  ) %>%
  filter(gSSMU %in% c("1", "2"), season %in% c("S", "W")) %>%
  group_by(cal.yr, gSSMU, season) %>%
  summarise(catch = sum(catch, na.rm = TRUE), .groups = "drop") %>%
  mutate(catch = replace_na(catch, 0))   # <-- NA to 0  # there is no missing data, but rather no fishing

## We impute data, based on our analysis for Median LKB -/+ SAM
source("./javier_analysis/survey_imputed2.R") # estimates 'survey' after filtering poor coverage


# --- Combine & compute fraction ---
catch_fraction <- inner_join(survey, fishery_by, by = c("cal.yr", "gSSMU", "season")) %>%
  mutate(
    fraction = if_else(!is.na(survey) & survey > 0,
                       round(catch / survey, 5), NA_real_)
  ) %>%
  arrange(gSSMU, cal.yr, season)


catch_imputed <- inner_join(LKB_imputed, fishery_by, by = c("cal.yr", "gSSMU", "season")) %>%
  mutate(
    fraction = if_else(!is.na(LKB) & LKB > 0,
                       round(catch / LKB, 5), NA_real_)
  ) %>%
  arrange(gSSMU, cal.yr, season)


# Total Krill Biomass surveyed per year
TotSurvey <- survey %>%
  filter(season == "S", gSSMU <3) %>%
  group_by(cal.yr) %>%
  summarise(Krill_biomass = sum(survey, na.rm = TRUE), .groups = "drop")
min(TotSurvey$Krill_biomass); max(TotSurvey$Krill_biomass)

# Save
write_csv(TotSurvey, file.path(out_dir, "Total_Krill_Survey_gSMMU_1_2.csv"))
write_csv(catch_fraction, file.path(out_dir, "LHR_fraction.csv"))
write_csv(catch_imputed, file.path(out_dir, "LHR_imputed.csv"))



### -------- PLOTTING ---------
# SURVEY DATA
# facet labels and ordering (put gSSMU2 on top)
facet_labs <- c(`1` = "Bransfield", `2` = "SSIW + EI")
facet_order <- c("2", "1")  # top to bottom in facet_wrap(ncol = 1)

p <- survey %>%
  filter(gSSMU %in% c("1", "2")) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)        # 2 on top, 1 below
  ) %>%
  ggplot(aes(x = cal.yr, y = survey, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue"), breaks = c("S", "W")) +
  facet_wrap(~ gSSMU, ncol = 1, labeller = as_labeller(facet_labs)) +  # shared y by default
  geom_point(data = LKB_imputed,
             aes(x = cal.yr, y = LKB, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  scale_y_continuous(labels = label_comma()) +  # e.g., 1,000,000 (no scientific)
  labs(
    title = "Krill Survey Biomass by gSSMU",
    x = "Year",
    y = "Biomass (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Krill_Survey.png")
ggsave(filename = outfile, plot = p, width = 10, height = 10, dpi = 300, bg = "white")

# Separate figures for gSSMU1 and 2
p <- survey %>%
  filter(gSSMU %in% c("1")) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)        # 2 on top, 1 below
  ) %>%
  ggplot(aes(x = cal.yr, y = survey, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue"), breaks = c("S", "W")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +  # shared y by default
  geom_point(data = filter(LKB_imputed, gSSMU == "1"),
             aes(x = cal.yr, y = LKB, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  scale_y_continuous(labels = label_comma()) +  # e.g., 1,000,000 (no scientific)
  labs(
    title = "Krill Survey - BS",
    x = "Year",
    y = "Biomass (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Krill_Survey_gSSMU1.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")


p <- survey %>%
  filter(gSSMU %in% c("2")) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)        # 2 on top, 1 below
  ) %>%
  ggplot(aes(x = cal.yr, y = survey, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue"), breaks = c("S", "W")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +  # shared y by default
  geom_point(data = filter(LKB_imputed, gSSMU == "2"),
             aes(x = cal.yr, y = LKB, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  scale_y_continuous(labels = label_comma()) +  # e.g., 1,000,000 (no scientific)
  labs(
    title = "Krill Survey - SSIW + EI",
    x = "Year",
    y = "Biomass (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Krill_Survey_gSMMU2.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")


# FISHERY DATA
# facet labels and ordering
facet_labs <- c(`1` = "Bransfield",
                `2` = "SSIW + EI")  
facet_order <- c("2", "1")  # SSIW + EI on top

p <- fishery_by %>%
  filter(gSSMU %in% c(1, 2)) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)
  ) %>%
  ggplot(aes(x = cal.yr, y = catch, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1, labeller = as_labeller(facet_labs)) +
  scale_y_continuous(labels = label_comma()) +   # show full digits (1,000,000)
  labs(
    title = "Krill Fishery Catches by gSSMU",
    x = "Year",
    y = "Catch (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "fishery_catch.png")
ggsave(filename = outfile, plot = p, width = 10, height = 10, dpi = 300, bg = "white")

# Figures separated by gSSMU
p <- fishery_by %>%
  filter(gSSMU %in% c(1)) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)
  ) %>%
  ggplot(aes(x = cal.yr, y = catch, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +
  scale_y_continuous(labels = label_comma()) +   # show full digits (1,000,000)
  labs(
    title = "Krill Fishery Catches - BS",
    x = "Year",
    y = "Catch (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "fishery_catch_gSSMU1.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")

p <- fishery_by %>%
  filter(gSSMU %in% c(2)) %>%
  mutate(
    gSSMU  = factor(gSSMU, levels = facet_order)
  ) %>%
  ggplot(aes(x = cal.yr, y = catch, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +
  scale_y_continuous(labels = label_comma()) +   # show full digits (1,000,000)
  labs(
    title = "Krill Fishery Catches - SSIW + EI",
    x = "Year",
    y = "Catch (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "fishery_catch_gSSMU2.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")


##### PLOT CATCH RATIO  ####

# facet labels
facet_labs <- c(`1` = "Bransfield", `2` = "SSIW + EI")

p <- catch_fraction %>%
  filter(gSSMU %in% c("1", "2")) %>%
  ggplot(aes(x = cal.yr, y = fraction, color = season, group = season)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1, labeller = as_labeller(facet_labs)) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),   # show as % with decimals
    limits = c(0, 0.32)                          # 0–100% range
  ) +
  geom_point(data = catch_imputed,
             aes(x = cal.yr, y = fraction, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  labs(
    title = "Fishery Catch as Fraction of Survey Biomass",
    x = "Year",
    y = "Catch / Biomass",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Catch_Ratio.png")
ggsave(filename = outfile, plot = p, width = 10, height = 10, dpi = 300, bg = "white")

# Figures separate for gSSMU
p <- catch_fraction %>%
  filter(gSSMU %in% c("1")) %>%
  ggplot(aes(x = cal.yr, y = fraction, color = season, group = season)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +
  geom_point(data = filter(catch_imputed, gSSMU == "1"),
             aes(x = cal.yr, y = fraction, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),   # show as % with decimals
    limits = c(0, 0.32)                          # 0–100% range
  ) +
  labs(
    title = "Local Harvest Ratio - BS",
    x = "Year",
    y = "Catch / Biomass",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Harvest_Ratio_gSSMU1.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")

p <- catch_fraction %>%
  filter(gSSMU %in% c("2")) %>%
  ggplot(aes(x = cal.yr, y = fraction, color = season, group = season)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_point(shape = 20, size = 5) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) + #, labeller = as_labeller(facet_labs)) +
  geom_point(data = filter(catch_imputed, gSSMU == "2"),
             aes(x = cal.yr, y = fraction, col = season, group = season),
             shape = 21,          # open circle
             size  = 3.5,
             na.rm = TRUE) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.05),   # show as % with decimals
    limits = c(0, 0.10)                          # 0–100% range
  ) +
  labs(
    title = "Harvest Ratio - SSIW + EI",
    x = "Year",
    y = "Catch / Biomass",
    color = "Season"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ); p
outfile <- file.path(out_dir, "Harvest_Ratio_gSSMU2.png")
ggsave(filename = outfile, plot = p, width = 6, height = 7.5, dpi = 300, bg = "white")



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

## End of Watters Code ##


