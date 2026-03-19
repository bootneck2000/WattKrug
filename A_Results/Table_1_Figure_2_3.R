### Watter et al. 2020

# Model imputed LKB
# Re-creation of equations (1) and (2): mean(k) and LKBij

library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(scales)
library(stringr)
library(patchwork)
library(knitr)
library(kableExtra)
#### 1. Import data from the model ####

# krill survey biomass
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)

# Filter low transect coverage (<10th percentil)
survey_filtered <- survey %>%
  filter(gSSMU %in% c(1, 2))
percentiles <- survey_filtered %>%
  group_by(gSSMU) %>%
  summarise(p10 = quantile(nmi.count, probs = 0.10, na.rm = TRUE))
percentiles$p10 <- round(percentiles$p10,0)   #  ceiling(percentiles$p10/10)*10
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

# Import ONI data
oni<-read.csv("./Supplementary Files/oni.csv",stringsAsFactors = FALSE)
oni$yr<-ifelse(is.element(oni$SEAS,c("OND","NDJ")),oni$YR+1,oni$YR)
oni$season<-ifelse(is.element(oni$SEAS,c("OND","NDJ","DJF","JFM")),"S",NA)
oni$season<-ifelse(is.element(oni$SEAS,c("AMJ","MJJ","JJA","JAS")),"W",oni$season)
oni<-na.omit(oni)
oni<-tapply(oni$ANOM,list(oni$yr,oni$season),mean)
oni<-data.frame(yr=rep(dimnames(oni)[[1]],2),season=rep(dimnames(oni)[[2]],each=dim(oni)[1]),oni=c(oni))
oni <- oni %>% mutate(cal.yr = as.integer(yr)) %>% dplyr::select(-yr)
oni$oni.class<-ifelse(oni$oni <= -0.5, "Cool","Neutral")
oni$oni.class<-ifelse(oni$oni >=0.5, "Warm",oni$oni.class)
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

fishery <- fishery %>% mutate(cal.yr = as.integer(cal.yr)) 
fishery$catch[is.na(fishery$catch)] <- 0

catches <- fishery %>%
  filter(gSSMU <3) %>%
  group_by(cal.yr) %>%
  summarise(Tot.catch = sum(catch))


# junk <- readRDS("./javier_analysis/Watters_model_output1/junk.rds")
# data <- junk %>% filter(season == "S") %>% filter(!is.na(survey)) %>% distinct(survey)

#### 2. Imputing LKB data based on SAM ####

# join variables
data <- left_join(fishery, survey, by = c("cal.yr", "gSSMU", "season", "matchme")) 
data <- data %>% filter(gSSMU %in% c(1,2))
data <- left_join(data, sam, by = c("cal.yr", "season")) 
data <- data %>% mutate(sam.sign = (ifelse(sam<0, "Neg", "Pos")))
data <- left_join(data, oni, by = c("cal.yr", "season")) 

## Estimating mean(logLKB) for survey data
data$LnSurvey <- log(data$survey)

# Survey data
#LKB_survey <- data 
Table.01 <- data %>%
  filter(!is.na(LnSurvey),
         gSSMU %in% c(1, 2),
         sam.sign %in% c("Neg", "Pos"),
         cal.yr < 2012) %>%
  group_by(gSSMU,sam.sign) %>%
  summarise(
    n       = n(),
    mu  = mean(LnSurvey),
    sd  = sd(LnSurvey),
    .groups = "drop"
  )
knitr::kable(Table.01, digits = 3)


# Imputing missing data
data <- data %>%
  mutate(
    # 0 if survey available, 1 if missing
    # impute.me = 1 # if_else(is.na(survey), 1, 0),  #Impute the whole series, for comparison
    
    # fill log survey
    LogSurvey_imputed = case_when(
      gSSMU == 1  & sam.sign == "Neg" ~ as.numeric(Table.01[1, 4]),  # remove from begining: 'impute.me == 1 &' 
      gSSMU == 1  & sam.sign == "Pos" ~ as.numeric(Table.01[2, 4]),
      gSSMU == 2  & sam.sign == "Neg" ~ as.numeric(Table.01[3, 4]),
      gSSMU == 2  & sam.sign == "Pos" ~ as.numeric(Table.01[4, 4]),
      TRUE ~ NA_real_
    ),
    
    # back-transform
    survey_imputed = exp(LogSurvey_imputed)
  )

#### 3. Statistics ####

### Statistical Analysis 
library(coin)

# Difference LKB between -/+SAM for gSSMU 1
g1 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 1) %>% filter(cal.yr<2012)
g1 <- g1 %>% mutate(sam.sign = as.factor(sam.sign))
g1.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g1, distribution = "exact")

# Difference LKB between -/+SAM for gSSMU 2
g2 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 2) %>% filter(cal.yr<2012)
g2 <- g2 %>% mutate(sam.sign = as.factor(sam.sign))
g2.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g2, distribution = "exact")

# Difference LKB between -/+SAM between gSSMUs
gNeg <- data %>% filter(!is.na(survey)) %>% 
  filter(sam.sign == "Neg") %>% filter(cal.yr<2012)
gNeg <- gNeg %>% mutate(gSSMU = as.factor(gSSMU))
gNeg.stat <- coin::oneway_test(LnSurvey ~ gSSMU, data = gNeg, distribution = "exact")

gPos <- data %>% filter(!is.na(survey)) %>% 
  filter(sam.sign == "Pos") %>% filter(cal.yr<2012)
gPos <- gPos %>% mutate(gSSMU = as.factor(gSSMU))
gPos.stat <- coin::oneway_test(LnSurvey ~ gSSMU, data = gPos, distribution = "exact")


#### Table 1 

Table.01$'p-value_1' <- NA
Table.01$`p-value_1`[1] <- pvalue(g1.stat)
Table.01$`p-value_1`[3] <- pvalue(g2.stat)

Table.01$Median_LKB <- exp(Table.01$mu)

Table.01$'p-value_2' <- NA
Table.01$`p-value_2`[1] <- pvalue(gNeg.stat)
Table.01$`p-value_2`[4] <- pvalue(gPos.stat)


knitr::kable(Table.01, digits = 3)

script_path <- "/Users/god/Documents/R workspace/WattKrug/A_Results/Table_1_Figure_2_3.R"

# Output directory = folder containing that script
out_dir <- "./A_Results/"

tab <- knitr::kable(Table.01, digits = 3, format = "html") |>
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "condensed"),
    full_width = FALSE,
    font_size = 14
  )

# Best: save to PDF (vector, publication-grade)
kableExtra::save_kable(
  tab,
  file = file.path(out_dir, "Table_1.pdf")
)
# Or save to PNG (raster)
kableExtra::save_kable(
  tab,
  file = file.path(out_dir, "Table_1.png"),
  density = 300
)

### 4. Plots to see relation between SAM and Survey data ####

ENV.LKB <- data %>%
  filter(!is.na(survey), cal.yr < 2012) %>%
  mutate(
    gSSMU = as.factor(gSSMU),
  ) %>%
  filter(season == "S") %>%
  dplyr::select(1:8)

LKB.imputed <- data %>%
    filter(cal.yr < 1996 | cal.yr > 2011) %>%
  mutate(gSSMU = as.factor(gSSMU)) %>%
  filter(season == "S") %>%
  dplyr::select(1:5, 7:8, 10)
LKB.imputed$LogSurvey_imputed = 0

Figure_2a <- ggplot(ENV.LKB, aes(x = sam, y = survey, col = gSSMU, group = gSSMU)) +
  # Reference lines
  geom_hline(yintercept = 1e+06, color = "black",   # 1Mt
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5) +  # SAM = 0
  geom_hline(yintercept = Table.01$Median_LKB[1], color = "red",      # mean LKB, gSSMU1, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Median_LKB[2], color = "red",      # mean LKB, gSSMU1, SAM Pos
             linewidth = 1) +
  geom_hline(yintercept = Table.01$Median_LKB[3], color = "blue",      # mean LKB, gSSMU2, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Median_LKB[4], color = "blue",      # mean LKB, gSSMU2, SAM Pos
             linewidth = 1) +
  # Points
  geom_point(alpha = 0.75, size = 3.5, na.rm = TRUE) +
  scale_color_manual(values = c("1" = "red", "2" = "blue"), 
                     breaks = c("1", "2"), name = "gSSMU") +
  # Labels and theme
  # y-axis: always show a tick at 1e6
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = union(pretty(ENV.LKB$survey, n = 6), 1e6) %>% sort()
  ) +
  labs(title = "A) LKB vs. SAM value",
       x = "SAM index", y = "Krill Biomass (ton)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
    ); Figure_2a

ENV.LKB <-  ENV.LKB %>% 
  mutate(
    survey   = as.numeric(survey),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU, levels = c(1, 2)),
  )

# --- Plot: bbplot per SAM sign, faceted by gSSMU ---
Figure_2b <- ggplot(ENV.LKB,
                    aes(x = sam.sign, y = survey, color = gSSMU)) +
  geom_hline(yintercept = 1e6, color = "black",                 # 1 Mt 
             linewidth = 1.1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Median_LKB[1], color = "red",  # mean LKB, gSSMU1, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Median_LKB[2], color = "red",  # mean LKB, gSSMU1, SAM Pos
             linewidth = 1) +
  geom_hline(yintercept = Table.01$Median_LKB[3], color = "blue", # mean LKB, gSSMU2, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Median_LKB[4], color = "blue", # mean LKB, gSSMU2, SAM Pos
             linewidth = 1) +
    geom_boxplot(width = 0.7,
               alpha = 0.6,
               outlier.shape = NA) +
  geom_jitter(width = 0.15,
              alpha = 0.9,
              size = 3,
              na.rm = TRUE) +
  facet_wrap(~ gSSMU, ncol = 2) +
  scale_color_manual(values = c("1" = "red", "2" = "blue"),
                     name = "gSSMU") +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = union(pretty(ENV.LKB$survey, n = 6), 1e6) |> sort()
  ) +
  labs(
    title = "B) LKB vs SAM sign",
    x = "SAM sign",
    y = "Krill Biomass (ton)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ); Figure_2b

# Save a stack figure
Fig_2_stack <- Figure_2a / Figure_2b

ggsave(
  filename = "./A_Results/Figure_2_stacked.png",
  plot     = Fig_2_stack,
  width    = 8,
  height   = 10,
  dpi      = 300
)

### 5. Plots: Survey data ####

### Figure 3: survey LKB vs. imputed data gSSMU1 and 2
Fig3a <- data %>%
  filter(gSSMU == "1") %>%
  ggplot(aes(x = cal.yr, y = survey/1000000, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(data = filter(survey, gSSMU == "1"),
             shape = 16, size = 4) +   # filled circles
  geom_point(data = filter(survey, gSSMU == "1", season == "W"),
             shape = 16, size = 4) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) +
  geom_point(data = filter(data, gSSMU == "1"),
             aes(x = cal.yr, y = survey_imputed/1000000, color = season, group = season),
             shape = 21, size = 3, na.rm = TRUE) +
  scale_y_continuous(labels = scales::label_comma(),
                     breaks = seq(0, max(data$survey/1000000, na.rm = TRUE), by = 1)) +
  labs(title = "Survey - BS",
       x = "Year", y = "Biomass (Million tons)", color = "Season") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig3a


Fig3b <- data %>%
  filter(gSSMU == "2") %>%
  ggplot(aes(x = cal.yr, y = survey/1000000, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(data = filter(survey, gSSMU == "2", cal.yr <= 2012),
             shape = 16, size = 4) +   # filled circles
  geom_point(data = filter(survey, gSSMU == "2", season == "W"),
             shape = 16, size = 4) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  facet_wrap(~ gSSMU, ncol = 1) +
  geom_point(data = filter(data, gSSMU == "2"),
             aes(x = cal.yr, y = survey_imputed/1000000, color = season, group = season),
             shape = 21, size = 3, na.rm = TRUE) +
  scale_y_continuous(labels = scales::label_comma(),
                     breaks = seq(0, max(data$survey/1000000, na.rm = TRUE), by = 1)) +
  labs(title = "Survey - SSWI+EI",
       x = "Year", y = "Biomass (Million tons)", color = "Season") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig3b


### 6. Plots: Fishery data ####

### Figure 4: Fishery vs. Year by gSSMU

# Figures separated by gSSMU
Fig4a <- data %>%
  filter(gSSMU == 1) %>%
  ggplot(aes(x = cal.yr, y = catch, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 16, size = 4) +
  facet_wrap(~ gSSMU, ncol = 1) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  scale_y_continuous(
    labels = scales::label_comma(),
    breaks = seq(0, max(data$catch, na.rm = TRUE), by = 1e4)
  ) +
  labs(
    title = "Fishing - BS",
    x = "Year",
    y = "Catch (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig4a

Fig4b <- data %>%
  filter(gSSMU == 2) %>%
  ggplot(aes(x = cal.yr, y = catch, color = season, group = season)) +
  geom_line(alpha = 0.7) +
  geom_point(shape = 16, size = 4) +
  facet_wrap(~ gSSMU, ncol = 1) +
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  scale_y_continuous(
    labels = scales::label_comma(),
    breaks = seq(0, max(data$catch, na.rm = TRUE), by = 1e4)
  ) +
  labs(
    title = "Fishing - SSWI+EI",
    x = "Year",
    y = "Catch (tons)",
    color = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig4b


### 7. LHR: estimations and plots ####

data <- data %>% mutate(
  LHR_obs = catch/survey,
  LHR_imp = catch/survey_imputed
  ) %>%
  dplyr::select(cal.yr, season, gSSMU, survey, catch, LHR_obs, sam, sam.sign,
         oni, oni.class, LnSurvey, LogSurvey_imputed, survey_imputed,
         LHR_imp)
  

### Figure 5: LHR

# Figures separated by gSSMU
Fig5a <- data %>%
  filter(gSSMU == 1) %>%
  ggplot(aes(x = cal.yr, group = season, color = season)) +
  
  # observed LHR
  geom_line(aes(y = LHR_obs), alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = LHR_obs), shape = 16, size = 3.5, na.rm = TRUE) +

  # LHR levels
  geom_hline(yintercept = 0.01, color = "black",  # very precautionary
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.10, color = "black",  # precautionary
             linewidth = 1) +
  
  # imputed LHR (dashed line, hollow points)
  geom_line(aes(y = LHR_imp), linetype = "dashed", alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = LHR_imp), shape = 21, 
             size = 3, stroke = 0.7, na.rm = TRUE) +
  
  facet_grid(season ~ gSSMU, scales = "free_y") +
  
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  
  scale_y_continuous(
    labels = label_comma()
  ) +
  
  labs(
    title = "LHR - BS",
    x = "Year",
    y = "Local Harvest Rate",
    color = "Season"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig5a

Fig5b <- data %>%
  filter(gSSMU == 2) %>%
  ggplot(aes(x = cal.yr, group = season, color = season)) +
  
  # observed LHR
  geom_line(aes(y = LHR_obs), alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = LHR_obs), shape = 16, size = 3.5, na.rm = TRUE) +
  
  # LHR levels
  geom_hline(yintercept = 0.01, color = "black",  # very precautionary
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = 0.10, color = "black",  # precautionary
             linewidth = 1) +
  
  # imputed LHR (dashed line, hollow points)
  geom_line(aes(y = LHR_imp), linetype = "dashed", alpha = 0.7, na.rm = TRUE) +
  geom_point(aes(y = LHR_imp), shape = 21, 
             size = 3, stroke = 0.7, na.rm = TRUE) +
  
  # facet_grid(season ~ gSSMU, scales = "free_y") +
  
  scale_color_manual(values = c("S" = "red", "W" = "blue")) +
  
  scale_y_continuous(
    labels = label_comma()
  ) +
  
  labs(
    title = "LHR - SSIW+EI",
    x = "Year",
    y = "Local Harvest Rate",
    color = "Season"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  ); Fig5b

# Save a stack figure with Figs 3, 4, 5
Fig3_stack <- (Fig3b + Fig4b + Fig5b) /
  (Fig3a + Fig4a + Fig5a)

ggsave(
  filename = "./A_Results/Figures_3_5_stacked.png",
  plot     = Fig3_stack,
  width    = 10,
  height   = 11,
  dpi      = 300
)
