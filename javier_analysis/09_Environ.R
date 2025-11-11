## Environmental Variables

## Watters Code:

# SOUTHERN ANNULAR MODE
sam<-read.csv("./Supplementary Files/sam.csv")
names(sam)<-c("yr","mo","sam")
sam$season<-ifelse(is.element(sam$mo,c(10:12,1:3)),"S","W")
sam$cal.yr<-ifelse(is.element(sam$mo,10:12),sam$yr+1,sam$yr)
sam<-tapply(sam$sam,list(sam$cal.yr,sam$season),mean)
sam<-data.frame(cal.yr=rep(dimnames(sam)[[1]],2),season=rep(dimnames(sam)[[2]],each=dim(sam)[1]),sam=c(sam))
# out$sam<-sam$sam[match(paste(out$cal.yr,out$season,sep="|"),paste(sam$YEAR,sam$season,sep="|"))]
# out$sam.sign<-ifelse(out$sam<0,"Neg","Pos")


# OCEANIC NINO INDEX
oni<-read.csv("./Supplementary Files/oni.csv",stringsAsFactors = FALSE)
oni$cal.yr<-ifelse(is.element(oni$SEAS,c("OND","NDJ")),oni$YR+1,oni$YR)
oni$season<-ifelse(is.element(oni$SEAS,c("OND","NDJ","DJF","JFM")),"S",NA)
oni$season<-ifelse(is.element(oni$SEAS,c("AMJ","MJJ","JJA","JAS")),"W",oni$season)
oni<-na.omit(oni)
oni<-tapply(oni$ANOM,list(oni$cal.yr,oni$season),mean)
oni<-data.frame(cal.yr=rep(dimnames(oni)[[1]],2),season=rep(dimnames(oni)[[2]],each=dim(oni)[1]),oni=c(oni))
# out$oni<-oni$oni[match(paste(out$cal.yr,out$season,sep="|"),paste(oni$yr,oni$season,sep="|"))]
# out$oni.class<-ifelse(out$oni <= -0.5, "Cool","Neutral")
# out$oni.class<-ifelse(out$oni >=0.5, "Warm",out$oni.class)
#
#


# ==================== ANALYSIS  ====================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})
suppressWarnings(RNGkind("default"))  # silence harmless RNG warnings on some R builds

# -------------------- Preconditions & setup ------------------------------
# source("./javier_analysis/Watters_Function.R") # run Watters' script
# junk <- make.localhr.data()  # CHANGE FOR SURVEY DATA FILTERED BY COVERAGE

source("./javier_analysis/survey_filtered.R") # estimates 'survey' after filtering poor coverage


# --- Output locations ---------------------------
out_dir <- "./javier_analysis/Env_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Add ENV data to survey:
oni <- oni %>% mutate(cal.yr =  as.integer(cal.yr))
sam <- sam %>% mutate(cal.yr =  as.integer(cal.yr))

ENV.LKB <- survey %>%
  left_join(oni, by = c("cal.yr", "season"), relationship = "many-to-one") %>%   # adds ONI columns
  left_join(sam, by = c("cal.yr", "season"), relationship = "many-to-one")       # adds SAM columns

# Extract imputed data from Supplementary Fig 7
dataS7 <- readRDS("./javier_analysis/Data_Fig_S7.rds")
S7.imputed <- dataS7 %>% filter(impute.me == 1) %>%
 select(cal.yr, gSSMU, survey, catch, sam, sam.sign, oni, oni.class, 
        impute.me) %>%
 distinct() 


# Make graphs to check relation between ONI/SAM and Survey data
ENV.LKB <- ENV.LKB %>%
  mutate(
    gSSMU     = as.factor(gSSMU),
    sam       = suppressWarnings(as.numeric(sam)),
    oni       = suppressWarnings(as.numeric(oni)),
    sam.sign  = as.factor(ifelse(sam<0, "Neg", "Pos")),
  ) %>%
  filter(season == "S")

S7.imputed <- S7.imputed %>%
  mutate(
    cal.yr    = as.integer(cal.yr),
    gSSMU     = as.factor(gSSMU),
    sam       = suppressWarnings(as.numeric(sam)),
    oni       = suppressWarnings(as.numeric(oni)),
    sam.sign  = as.factor(sam.sign),
    oni.class = as.factor(oni.class)
  )

# (i) x = sam, y = survey
fig1a <- ggplot(ENV.LKB, aes(x = sam, y = survey, col = gSSMU, group = gSSMU)) +
  # Reference lines
  geom_hline(yintercept = 1e+06, color = "grey50", linewidth = 1.5) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5) +
  # Points
  geom_point(alpha = 0.75, size = 2.5, na.rm = TRUE) +
  # Add imputed data
  geom_point(
    data = S7.imputed,
    aes(x = sam, y = survey, col = gSSMU, group = gSSMU),
    shape = 2,          # filled circle (optional)
    size  = 3.5,
    na.rm = TRUE
  ) +
  # Facet and colors
  # facet_wrap(~ gSSMU, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("1" = "red", "2" = "blue"), 
                     breaks = c("1", "2"), name = "gSSMU") +
  # Labels and theme
  labs(title = "survey ~ sam (faceted by gSSMU)", x = "sam", y = "survey") +
  theme_minimal(base_size = 14); fig1a

outfile <- file.path(out_dir, "Survey_vs_SAM.png")
ggsave(filename = outfile, plot = fig1a, width = 10, height = 10, dpi = 300, bg = "white")


# (ii) x = sam.sign, y = survey
fig2a <- ggplot(ENV.LKB, aes(x = sam.sign, y = survey, col = gSSMU, group = gSSMU)) +
  # Reference lines
  geom_hline(yintercept = 1e+06, color = "grey50", linewidth = 1.5) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5) +
  # Points
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, na.rm = TRUE) +
  # Add imputed data
  geom_jitter(
    data = S7.imputed,
    aes(x = sam.sign, y = survey, col = gSSMU, group = gSSMU),
    width = 0.2,
    shape = 2,          # filled circle (optional)
    size  = 3.5,
    na.rm = TRUE
  ) +
  # Facet and colors
  # facet_wrap(~ gSSMU, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("1" = "red", "2" = "blue"), 
                     breaks = c("1", "2"), name = "gSSMU") +
  # stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "black") +
  # Labels and theme
  scale_y_continuous(
    labels = scales::label_number(scale = 1, big.mark = ","),
    breaks = pretty(ENV.LKB$survey, n = 6) %>% union(1e6) %>% sort()
  ) +
  labs(title = "survey ~ sam.sign (faceted by gSSMU)", x = "sam.sign", y = "survey") +
  theme_minimal(base_size = 14); fig2a

outfile <- file.path(out_dir, "Survey_vs_SAM_sign.png")
ggsave(filename = outfile, plot = fig2a, width = 7, height = 7, dpi = 300, bg = "white")


# (iii) Boxplot: x = sam.sign + gSSMU, y = survey biomass

# --- Prep: bind rows with a dataset flag, keep types tidy ---
combined <- bind_rows(
  ENV.LKB     %>% mutate(source = "Original"),
  S7.imputed  %>% mutate(source = "Imputed")
) %>%
  mutate(
    survey   = as.numeric(survey),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU, levels = c(1, 2)),
    source   = factor(source, levels = c("Original", "Imputed"))
  )

# --- Plot: boxes per SAM sign, dodged by dataset, faceted by gSSMU ---
fig3 <- ggplot(combined, aes(x = sam.sign, y = survey,
                                fill = source, color = gSSMU)) +
  # reference line at 1,000,000
  geom_hline(yintercept = 1e6, color = "grey40", linewidth = 1.1, linetype = "dashed") +
  # boxplots: side-by-side Original vs Imputed within each SAM sign
  geom_boxplot(position = position_dodge(width = 0.85),
               width = 0.7, alpha = 0.7, outlier.shape = NA) +
  # (optional) raw points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85),
             alpha = 0.8, size = 3, na.rm = TRUE) +
  facet_wrap(~ gSSMU, ncol = 2) +
  # colors: gSSMU outlines red/blue; fills indicate dataset
  scale_color_manual(values = c("1" = "red", "2" = "blue"), name = "gSSMU") +
  scale_fill_manual(values = c("Original" = "white", "Imputed" = "gray"),
                    name = "Dataset") +
  
  # y-axis: always show a tick at 1e6
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = union(pretty(combined$survey, n = 6), 1e6) %>% sort()
  ) +
  labs(title = "Survey data vs. SAM vs Imputed",
       x = "SAM sign", y = "Krill Biomass (ton)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ); fig3

outfile <- file.path(out_dir, "Survey_Imputed_vs_SAM_sign.png")
ggsave(filename = outfile, plot = fig3, width = 7, height = 7, dpi = 300, bg = "white")


#### Statistical Analysis
g1 <- ENV.LKB %>% filter(gSSMU == 1) # %>% filter(survey > 0)
coin::oneway_test(survey ~ sam.sign, data = g1, distribution = "approximate")
# Z = -1.1643, p-value = 0.2639
g1_summary <- g1 %>%
  group_by(sam.sign) %>%
  summarise(
    mean_survey   = mean(survey, na.rm = TRUE),
    median_survey = median(survey, na.rm = TRUE),
    n             = n()
  )
g1_summary

g2 <- ENV.LKB %>% filter(gSSMU == 2) # %>% filter(survey > 0)
coin::oneway_test(survey ~ sam.sign, data = g2, distribution = "approximate")
# Z = -1.9038, p-value = 0.0527
g2_summary <- g2 %>%
  group_by(sam.sign) %>%
  summarise(
    mean_survey   = mean(survey, na.rm = TRUE),
    median_survey = median(survey, na.rm = TRUE),
    n             = n()
  )
g2_summary
