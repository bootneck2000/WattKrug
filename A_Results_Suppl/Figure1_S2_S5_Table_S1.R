### Watter et al. 2020

# Model imputed LKB
# Re-creation of equations (1) and (2): mean(k) and LKBij

library(dplyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(scales)
library(stringr)

#### 1. Import data from the model ####

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

#### 2. Imputing LKB data based on SAM ####

# join variables
data <- left_join(fishery, survey, by = c("cal.yr", "gSSMU", "season", "matchme")) 
data <- data %>% filter(gSSMU %in% c(1,2))
data <- left_join(data, sam, by = c("cal.yr", "season")) 
data <- data %>% mutate(sam.sign = (ifelse(sam<0, "Neg", "Pos")))

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
      gSSMU == 1  & sam.sign == "Neg" ~ as.numeric(Table.01[1, 4]),  # remove from beggining: 'impute.me == 1 &' 
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

g1 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 1) %>% filter(cal.yr<2012)
g1 <- g1 %>% mutate(sam.sign = as.factor(sam.sign))
g1.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g1, distribution = "approximate")

g2 <- data %>% filter(!is.na(survey)) %>% 
  filter(gSSMU == 2) %>% filter(cal.yr<2012)
g2 <- g2 %>% mutate(sam.sign = as.factor(sam.sign))
g2.stat <- coin::oneway_test(LnSurvey ~ sam.sign, data = g2, distribution = "approximate")

#### Table 1 

Table.01$'p-value' <- NA
Table.01$`p-value`[1] <- pvalue(g1.stat)
Table.01$`p-value`[3] <- pvalue(g2.stat)

Table.01$Mean_LKB <- exp(Table.01$mu)

knitr::kable(Table.01, digits = 3)

### 4. Plots to see relation between SAM and Survey data ####

ENV.LKB <- data %>%
  filter(!is.na(survey), cal.yr < 2012) %>%
  mutate(
    gSSMU = as.factor(gSSMU),
  ) %>%
  filter(season == "S") %>%
  select(1:8)

LKB.imputed <- data %>%
    filter(cal.yr < 1996 | cal.yr > 2011) %>%
  mutate(gSSMU = as.factor(gSSMU)) %>%
  filter(season == "S") %>%
  select(1:5, 7:8, 10)
LKB.imputed$LogSurvey_imputed = 0
# LKB.imputed$survey_imputed = 0

# Estimate imputed values for each year, according to gSSMU and SAM sign
set.seed(123)  # optional, for reproducibility

# log limits, as established by Watters et al.
lower_log <- log(10000)
upper_log <- log(100000000)

for(i in 1:length(LKB.imputed$cal.yr)) {
  if(LKB.imputed$gSSMU[i] == "1") {
    if(LKB.imputed$sam.sign[i] == "Neg") {
      mu_ls <- Table.01$mu[1]
      sd_ls <- Table.01$sd[1]
    }
    if(LKB.imputed$sam.sign[i] == "Pos") {
      mu_ls <- Table.01$mu[2]
      sd_ls <- Table.01$sd[2]
    }
  }
  if(LKB.imputed$gSSMU[i] == "2") {
    if(LKB.imputed$sam.sign[i] == "Neg") {
      mu_ls <- Table.01$mu[3]
      sd_ls <- Table.01$sd[3]
    }
    if(LKB.imputed$sam.sign[i] == "Pos") {
      mu_ls <- Table.01$mu[4]
      sd_ls <- Table.01$sd[4]
    }
  }
  LKB.imputed$LogSurvey_imputed[i] <- rnorm(1, mu_ls, sd_ls)
  # truncate to [log(10000), log(1e8)]
  LKB.imputed$LogSurvey_imputed[i] <- pmin(pmax(LKB.imputed$LogSurvey_imputed[i], lower_log), upper_log)
  LKB.imputed$survey[i] <- exp(LKB.imputed$LogSurvey_imputed[i])
}


# (i) x = sam, y = survey
Figure_2a <- ggplot(ENV.LKB, aes(x = sam, y = survey, col = gSSMU, group = gSSMU)) +
  # Reference lines
  geom_hline(yintercept = 1e+06, color = "black",   # 1Mt
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 1.5) +  # SAM = 0
  geom_hline(yintercept = Table.01$Mean_LKB[1], color = "red",      # mean LKB, gSSMU1, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Mean_LKB[2], color = "red",      # mean LKB, gSSMU1, SAM Pos
             linewidth = 1) +
  geom_hline(yintercept = Table.01$Mean_LKB[3], color = "blue",      # mean LKB, gSSMU2, SAM Neg
             linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = Table.01$Mean_LKB[4], color = "blue",      # mean LKB, gSSMU2, SAM Pos
             linewidth = 1) +
  # Points
  geom_point(alpha = 0.75, size = 3.5, na.rm = TRUE) +
  # Add imputed data
  geom_point(
    data = LKB.imputed,
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
  # y-axis: always show a tick at 1e6
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = union(pretty(ENV.LKB$survey, n = 6), 1e6) %>% sort()
  ) +
  labs(title = "A) Survey data vs. SAM vs Imputed",
       x = "SAM index", y = "Krill Biomass (ton)") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
    ); Figure_2a

# --- Prep: bind rows with a dataset flag, keep types tidy ---
combined <- bind_rows(
  ENV.LKB     %>% mutate(source = "Original"),
  LKB.imputed  %>% mutate(source = "Imputed")
) %>%
  mutate(
    survey   = as.numeric(survey),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU, levels = c(1, 2)),
    source   = factor(source, levels = c("Original", "Imputed"))
  )

# --- Plot: boxes per SAM sign, dodged by dataset, faceted by gSSMU ---
Figure_2b <- ggplot(combined, 
                    aes(x = sam.sign, y = survey,
                        color = gSSMU)) +
  
  geom_hline(yintercept = 1e6, color = "black", linewidth = 1.1, linetype = "dashed") +
  
  # --- BOX PLOTS separated by source (Original/Imputed) ---
  geom_boxplot(aes(fill = source,
                   group = interaction(sam.sign, source)),
               position = position_dodge(width = 0.85),
               width = 0.7,
               alpha = 0.6,
               outlier.shape = NA) +
  
  # --- POINTS with separate shapes by source ---
  geom_jitter(aes(shape = source,
                  group  = interaction(sam.sign, source)),
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width = 0.85),
              alpha = 0.9,
              size = 3,
              na.rm = TRUE) +
  
  facet_wrap(~ gSSMU, ncol = 2) +
  
  # Colors for gSSMU (outline)
  scale_color_manual(values = c("1" = "red", "2" = "blue"),
                     name = "gSSMU") +
  
  # Fill for boxplots
  scale_fill_manual(values = c("Original" = "white", "Imputed" = "gray90"),
                    name = "Dataset") +
  
  # Shapes for points: filled circle vs open circle
  scale_shape_manual(values = c("Original" = 16,  # filled circle
                                "Imputed"  = 1),  # open circle
                     name = "Dataset") +
  
  scale_y_continuous(
    labels = scales::label_number(big.mark = ","),
    breaks = union(pretty(combined$survey, n = 6), 1e6) %>% sort()
  ) +
  
  labs(
    title = "B) Survey distribution by SAM sign, Dataset, and gSSMU",
    x = "SAM sign",
    y = "Krill Biomass (ton)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ); Figure_2b



### 5. Plots: Survey data ####

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




### 6. Estimating Frequency distributions for Krill Survey Data ####

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
curve_data <- Table.01 %>%
  mutate(
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU)
  ) %>%
  rowwise() %>%
  mutate(
    x = list(seq(mu - 4*sd, mu + 4*sd, length.out = 200))
  ) %>%
  unnest(x) %>%
  mutate(
    density = dnorm(x, mean = mu, sd = sd)
  ) %>%
  ungroup()


# 3. Plot: histogram of data + normal line from params_ln
x_min <- floor(min(hist_data$LnSurvey, na.rm = TRUE))
x_max <- ceiling(max(hist_data$LnSurvey, na.rm = TRUE))
x_ll <- seq(x_min, x_max, 0.5)
survey_hist <- ggplot(hist_data, aes(x = LnSurvey)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.5,
                 color = "black",
                 fill  = "grey80") +
  geom_line(data = curve_data,
            aes(x = x, y = density),
            linewidth = 1) +
  facet_grid(sam.sign ~ gSSMU) +
  scale_x_continuous(limits = c(x_min, x_max),
                     breaks = x_ll, 
                     labels = x_ll) +
  labs(x = "LnSurvey", y = "Density") +
  theme_minimal()
survey_hist

### 7. Estimating Frequency distributions for original IMPUTED LKB ####

### We estimate LKB imputed data using Equation (1) and (2) (Waters et al. 2020)

## Equation (2)
# K = U(0.1k, 10k); k = mean(log(LKB)), by SAM sign
# sigma = U(0.1s, 10s); s = sd(log(LKB)), by SAM sign

# First, creating mean(k) and sd(k) matrix
meanlogsummer <- matrix(
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("gSSMU1", "gSSMU2"),
    c("Neg", "Pos")
  )
)
meanlogsummer[1,] <- Table.01$mu[1:2] 
meanlogsummer[2,] <- Table.01$mu[3:4] 

sdlogsummer <- matrix(
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("gSSMU1", "gSSMU2"),
    c("Neg", "Pos")
  )
)
sdlogsummer[1,] <- Table.01$sd[1:2] 
sdlogsummer[2,] <- Table.01$sd[3:4] 

# Produce random imputed LKB values for each set o condition (gSSMU x SAM)
set.seed(123)  # optional, for reproducibility

n_sim <- 100000

# meanlogsummer and sdlogsummer are 2x2 matrices:
# rows: gSSMU1, gSSMU2
# cols: Neg, Pos

# log limits
lower_log <- log(10000)
upper_log <- log(100000000)

# prepare empty matrix for simulated LKB
sim_mat_log <- matrix(NA_real_, nrow = n_sim, ncol = 4)
colnames(sim_mat_log) <- c("g1.Neg", "g1.Pos", "g2.Neg", "g2.Pos")

# optional: store mu.ls and sd.ls if you want to keep them
mu_mat <- sim_mat_log
sd_mat <- sim_mat_log

col_index <- 1
for (i in 1:2) {          # gSSMU: 1,2
  for (j in 1:2) {        # sam.sign: Neg, Pos
    # draw mu.ls and sd.ls from uniform distributions
    mu_ls <- runif(
      n_sim,
      min = 0.1 * meanlogsummer[i, j],
      max = 10  * meanlogsummer[i, j]
    )
    sd_ls <- runif(
      n_sim,
      min = 0.1 * sdlogsummer[i, j],
      max = 10  * sdlogsummer[i, j]
    )
    
    # normal draws on log scale
    log_vals <- rnorm(
      n    = n_sim,
      mean = mu_ls,
      sd   = sd_ls
    )
    
    # truncate to [log(10000), log(1e8)]
    log_vals <- pmin(pmax(log_vals, lower_log), upper_log)
    
    # store
    sim_mat_log[, col_index] <- log_vals
    mu_mat[, col_index]      <- mu_ls
    sd_mat[, col_index]      <- sd_ls
    
    col_index <- col_index + 1
  }
}


# back to original (normal) scale
LKB.imp <- as.data.frame(exp(sim_mat_log))

str(LKB.imp)
head(LKB.imp)

# Plot histogram of imputed mean-log(LKB):

# Transform data to be used:
sim_long <- sim_mat_log %>%
  as.data.frame() %>%
  mutate(sim_id = row_number()) %>%
  pivot_longer(
    cols = c(g1.Neg, g1.Pos, g2.Neg, g2.Pos),
    names_to = c("gSSMU", "sam.sign"),
    names_sep = "\\."
  ) %>%
  rename(LnSurvey = value) %>%
  mutate(
    gSSMU    = factor(gSSMU, levels = c("g1", "g2"), labels = c(1, 2)),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
  )

#  1. Data for histograms
hist_imp_data <- sim_long %>%
  mutate(
    gSSMU    = factor(gSSMU, levels = c(1, 2)),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
  )

# 2. Data for fitted normal curves, built only from params_ln
# Imputed LKB 
Table.02 <- sim_long %>%
  filter(!is.na(LnSurvey),
         gSSMU %in% c(1, 2),
         sam.sign %in% c("Neg", "Pos")) %>%
  group_by(gSSMU,sam.sign) %>%
  summarise(
    n       = n(),
    mu  = mean(LnSurvey),
    sd  = sd(LnSurvey),
    .groups = "drop"
  )
knitr::kable(Table.02, digits = 3)

curve_data_imp <- Table.02 %>%
  mutate(
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU)
  ) %>%
  rowwise() %>%
  mutate(
    x = list(seq(mu - 4*sd, mu + 4*sd, length.out = 200))
  ) %>%
  unnest(x) %>%
  mutate(
    density = dnorm(x, mean = mu, sd = sd)
  ) %>%
  ungroup()

# 3. Plot: histogram of data + normal line from params_ln
x_min <- floor(min(hist_imp_data$LnSurvey, na.rm = TRUE))
x_max <- ceiling(max(hist_imp_data$LnSurvey, na.rm = TRUE))
x_ll <- seq(x_min, x_max, 1)
imputed_survey_hist <- ggplot(hist_imp_data, aes(x = LnSurvey)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 color = "black",
                 fill  = "grey80") +
  geom_line(data = curve_data_imp,
            aes(x = x, y = density),
            linewidth = 1) +
  facet_grid(sam.sign ~ gSSMU) +
  scale_x_continuous(limits = c(x_min, x_max),
                     breaks = x_ll, 
                     labels = x_ll) +
  labs(x = "LnSurvey", y = "Density") +
  theme_minimal()
imputed_survey_hist


### 8. Estimating Frequency distributions for imputed mean-log(LKB) - variation####

### We estimate LKB imputed data using only Equation (1)Waters et al. 2020)
# We use only mean-log(lkb) data; no Uniform distribution

# First, creating mean(k) and sd(k) matrix
meanlogsummer <- matrix(
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("gSSMU1", "gSSMU2"),
    c("Neg", "Pos")
  )
)
meanlogsummer[1,] <- Table.01$mu[1:2] 
meanlogsummer[2,] <- Table.01$mu[3:4] 

sdlogsummer <- matrix(
  nrow = 2,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("gSSMU1", "gSSMU2"),
    c("Neg", "Pos")
  )
)
sdlogsummer[1,] <- Table.01$sd[1:2] 
sdlogsummer[2,] <- Table.01$sd[3:4] 

# Produce random imputed LKB values for each set o condition (gSSMU x SAM)
set.seed(123)  # optional, for reproducibility

n_sim <- 1000

# meanlogsummer and sdlogsummer are 2x2 matrices:
# rows: gSSMU1, gSSMU2
# cols: Neg, Pos

# log limits
lower_log <- log(10000)
upper_log <- log(100000000)

# prepare empty matrix for simulated LKB
sim_mat_log2 <- matrix(NA_real_, nrow = n_sim, ncol = 4)
colnames(sim_mat_log2) <- c("g1.Neg", "g1.Pos", "g2.Neg", "g2.Pos")

# optional: store mu.ls and sd.ls if you want to keep them
mu_mat <- sim_mat_log
sd_mat <- sim_mat_log

col_index <- 1
for (i in 1:2) {          # gSSMU: 1,2
  for (j in 1:2) {        # sam.sign: Neg, Pos
    # draw mu.ls and sd.ls from meanlogsummer / sdlogsummer
    mu_ls <- meanlogsummer[i,j]
    sd_ls <- sdlogsummer[i,j]
    
    # normal draws on log scale
    log_vals <- rnorm(
      n    = n_sim,
      mean = mu_ls,
      sd   = sd_ls
    )
    
    # truncate to [log(10000), log(1e8)]
    log_vals <- pmin(pmax(log_vals, lower_log), upper_log)
    
    # store
    sim_mat_log2[, col_index] <- log_vals
    mu_mat[, col_index]      <- mu_ls
    sd_mat[, col_index]      <- sd_ls
    
    col_index <- col_index + 1
  }
}


# back to original (normal) scale
LKB.imp.mean <- as.data.frame(exp(sim_mat_log2))
head(LKB.imp.mean)

# Plot histogram of imputed mean-log(LKB):

# Transform data to be used:
sim_long2 <- sim_mat_log2 %>%
  as.data.frame() %>%
  mutate(sim_id = row_number()) %>%
  pivot_longer(
    cols = c(g1.Neg, g1.Pos, g2.Neg, g2.Pos),
    names_to = c("gSSMU", "sam.sign"),
    names_sep = "\\."
  ) %>%
  rename(LnSurvey = value) %>%
  mutate(
    gSSMU    = factor(gSSMU, levels = c("g1", "g2"), labels = c(1, 2)),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
  )

#  1. Data for histograms
hist_imp_data2 <- sim_long2 %>%
  mutate(
    gSSMU    = factor(gSSMU, levels = c(1, 2)),
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
  )

# 2. Data for fitted normal curves, built only from params_ln
# Imputed LKB 
Table.03 <- sim_long2 %>%
  filter(!is.na(LnSurvey),
         gSSMU %in% c(1, 2),
         sam.sign %in% c("Neg", "Pos")) %>%
  group_by(gSSMU,sam.sign) %>%
  summarise(
    n       = n(),
    mu  = mean(LnSurvey),
    sd  = sd(LnSurvey),
    .groups = "drop"
  )
knitr::kable(Table.03, digits = 3)

curve_data_imp2 <- Table.03 %>%
  mutate(
    sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
    gSSMU    = factor(gSSMU)
  ) %>%
  rowwise() %>%
  mutate(
    x = list(seq(mu - 4*sd, mu + 4*sd, length.out = 200))
  ) %>%
  unnest(x) %>%
  mutate(
    density = dnorm(x, mean = mu, sd = sd)
  ) %>%
  ungroup()

# 3. Plot: histogram of data + normal line from params_ln
x_min <- floor(min(hist_imp_data2$LnSurvey, na.rm = TRUE))
x_max <- ceiling(max(hist_imp_data2$LnSurvey, na.rm = TRUE))
x_ll <- seq(x_min, x_max, 1)
imputed_survey_hist.mean <- ggplot(hist_imp_data2, aes(x = LnSurvey)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 color = "black",
                 fill  = "grey80") +
  geom_line(data = curve_data_imp2,
            aes(x = x, y = density),
            linewidth = 1) +
  facet_grid(sam.sign ~ gSSMU) +
  scale_x_continuous(limits = c(x_min, x_max),
                     breaks = x_ll, 
                     labels = x_ll) +
  labs(x = "LnSurvey", y = "Density") +
  theme_minimal()
imputed_survey_hist.mean




