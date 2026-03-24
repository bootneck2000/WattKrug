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
percentiles$p10 <- round(percentiles$p10,0)
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

script_path <- "./A_Results_Suppl"

tab <- knitr::kable(Table.01, digits = 3, format = "html") |>
  kableExtra::kable_styling(
    bootstrap_options = c("striped", "condensed"),
    full_width = FALSE,
    font_size = 14
  )

# Best: save to PDF (vector, publication-grade)
kableExtra::save_kable(
  tab,
  file = file.path(script_path, "Table_S1.pdf")
)
# Or save to PNG (raster)
kableExtra::save_kable(
  tab,
  file = file.path(script_path, "Table_S1.png"),
  density = 300
)

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
hist_data2 <- hist_data %>% dplyr::select(gSSMU, sam.sign, LnSurvey)

# Testing for (log) normal distribution
shapiro_results1 <- hist_data2 %>%
  group_by(gSSMU, sam.sign) %>%
  summarise(
    n = n(),
    W       = shapiro.test(LnSurvey)$statistic,
    p_value = shapiro.test(LnSurvey)$p.value
  )

# Plot LKB frequecny, by gSSMU and SAM sign
# ---- p1: summer biomass ----
summer_obs <- hist_data %>%
    dplyr::select(-LogSurvey_imputed, -survey_imputed)

summer_obs$catch[is.na(summer_obs$catch)] <- 0
summer_obs$hr = summer_obs$catch/summer_obs$survey

p1_obs <- ggplot(summer_obs, aes(x = survey / 1e6,
                                     y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "coral3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_grid(sam.sign ~ gSSMU) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 0.15)) +
  scale_x_continuous(breaks = seq(0, 25, by = 1)) +
  labs(
    title    = "Observed summer krill biomass",
    x        = "Biomass (million tonnes)",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p1_obs

# ---- p1b: log biomass ----
p1b_obs <- ggplot(summer_obs, aes(x = log(survey),
                                      y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "#b2182b", color = "white", alpha = 0.85) +
  geom_vline(xintercept = log(1e6), linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_grid(sam.sign ~ gSSMU) +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(
    title    = "Observed summer krill biomass (log scale)",
    x        = "ln(Biomass in tonnes)",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p1b_obs

# ---- p2: summer harvest rate ----
p2_obs <- ggplot(summer_obs, aes(x = hr,
                                     y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 0.005, boundary = 0, fill = "coral3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_vline(xintercept = 0.10, linetype = "solid",  color = "black", linewidth = 0.7) +
  facet_grid(sam.sign ~ gSSMU) +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.30)) +
  scale_x_continuous(breaks = seq(0, 0.25, by = 0.01)) +
  labs(
    title    = "Observed summer harvest rate",
    x        = "Harvest Rate",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p2_obs

# Assigning LKB ad LHR categories
summer_cat <- summer_obs %>%
  mutate(
    bclass  = factor(ifelse(survey <= 1e6, "≤ 1 Mt", "> 1 Mt"),
                         levels = c("≤ 1 Mt", "> 1 Mt")),
    hrclass = factor(case_when(
      hr <= 0.01 ~ "≤ 0.01",
      hr >= 0.10 ~ "≥ 0.10",
      TRUE       ~ "0.01 – <0.10"),
      levels = c("≤ 0.01", "0.01 – <0.10", "≥ 0.10"))
  )
# ---- p3a observed: biomass class ----
p3a_obs <- summer_cat %>%
  count(gSSMU, sam.sign, bclass) %>%
  ggplot(aes(x = bclass, y = n / sum(n), fill = bclass)) +
  geom_col(width = 0.55, color = "white") +
  facet_grid(sam.sign ~ gSSMU) +
  coord_cartesian(ylim = c(0, 0.30)) +
  scale_fill_manual(values = c("≤ 1 Mt" = "coral3", "> 1 Mt" = "#b2182b")) +
  labs(
    title    = "Biomass class distribution — observed data",
    x        = "Biomass class",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none"); p3a_obs

# ---- p3b observed: harvest rate class ----
p3b_obs <- summer_cat %>%
  count(gSSMU, sam.sign, hrclass, .drop = FALSE) %>%
  ggplot(aes(x = hrclass, y = n / sum(n), fill = hrclass)) +
  geom_col(width = 0.55, color = "white") +
  facet_grid(sam.sign ~ gSSMU) +
  coord_cartesian(ylim = c(0, 0.30)) +
  scale_fill_manual(values = c(
    "≤ 0.01"       = "#f4a582",
    "0.01 – <0.10" = "#d6604d",
    "≥ 0.10"       = "#b2182b"
  )) +
  labs(
    title    = "Harvest rate class distribution — observed data",
    x        = "Harvest rate class",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none"); p3b_obs

out.dir <- "./A_Results_Suppl/Model01_filtered/"
ggsave(file.path(out.dir, "obs_freq_dist_summer_biomass.png"),     p1_obs,  width=8, height=6, dpi=300)
ggsave(file.path(out.dir, "obs_freq_dist_summer_biomass_log.png"), p1b_obs, width=8, height=6, dpi=300)
ggsave(file.path(out.dir, "obs_freq_dist_summer_hr.png"),          p2_obs,  width=8, height=6, dpi=300)
ggsave(file.path(out.dir, "obs_freq_dist_bclass.png"),             p3a_obs, width=8, height=6, dpi=300)
ggsave(file.path(out.dir, "obs_freq_dist_hrclass.png"),            p3b_obs, width=8, height=6, dpi=300)


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
x_min <- floor(min(hist_data2$LnSurvey, na.rm = TRUE))
x_max <- ceiling(max(hist_data2$LnSurvey, na.rm = TRUE))
x_ll <- seq(x_min, x_max, 0.5)
survey_hist <- ggplot(hist_data2, aes(x = LnSurvey)) +
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
survey_hist  # Figure S2





### 7. Estimating Frequency distributions for original IMPUTED LKB ####

out.dir <- "./A_Results/Watters_model/Imputed/"

hr.biomass.post <- readRDS(file.path(out.dir, "hr_biomass_post.rds"))
hr.biomass.mat <- readRDS(file.path(out.dir, "hr_biomass_mat.rds"))
hr.biomass.summ <- readRDS(file.path(out.dir, "hr_biomass_summ.rds"))
junk <- readRDS(file.path(out.dir, file = "junk.rds"))

# ============================================================
# Metadata tables (link MCMC node indices back to gSSMU & SAM)
# ============================================================

# Summer subset: indices 1:nsummerobs → rows of junk where season == "S"
summer_meta <- junk %>%
  filter(season == "S") %>%
  mutate(
    idx       = row_number(),
    impute.me = as.integer(is.na(survey)),
    gSSMU_lab = ifelse(gSSMU == 1, "gSSMU 1 (Bransfield)", "gSSMU 2 (Drake Passage)"),
    sam_lab   = ifelse(sam.sign == "Neg", "SAM Negative", "SAM Positive")
  )

# Full dataset: indices 1:nobs → all rows of junk (for bclass / hrclass)
all_meta <- junk %>%
  mutate(
    idx       = row_number(),
    impute.me = as.integer(is.na(survey) & season == "S"),
    gSSMU_lab = ifelse(gSSMU == 1, "gSSMU 1 (Bransfield)", "gSSMU 2 (Drake Passage)"),
    sam_lab   = ifelse(sam.sign == "Neg", "SAM Negative", "SAM Positive")
  )

# Imputed observation indices (summer only)
imputed_idx <- summer_meta %>% filter(impute.me == 1) %>% pull(idx)

# ============================================================
# Helper: extract node type and reshape to long format
# ============================================================
extract_long <- function(mat, node_type, meta, imputed_only = TRUE) {
  cols <- grep(paste0("^", node_type, "\\["), colnames(mat), value = TRUE)
  df <- mat[, cols, drop = FALSE] %>%
    as.data.frame() %>%
    mutate(sample = row_number()) %>%
    pivot_longer(cols = -sample, names_to = "node", values_to = "value") %>%
    mutate(idx = as.integer(str_extract(node, "\\d+"))) %>%
    left_join(meta %>% dplyr::select(idx, gSSMU_lab, sam_lab, impute.me), by = "idx")
  if (imputed_only) df <- filter(df, impute.me == 1)
  df
}

n_iter    <- nrow(hr.biomass.mat)   # 6,000 (2,000/chain × 3 chains)
n_imputed <- length(imputed_idx)

# ============================================================
# Plot 1: Imputed summer krill biomass
# ============================================================
p1_data <- extract_long(hr.biomass.mat, "summer", summer_meta, imputed_only = TRUE)


# 1. Data for histograms
hist_data_imp <- p1_data %>%
  mutate(
    sam.sign = factor(sam_lab),
    gSSMU    = factor(gSSMU_lab),
    LnSurvey = log(value)
  )
hist_data_imp <- hist_data_imp %>% select(gSSMU, sam.sign, LnSurvey)

# Testing for (log) normal distribution
library(nortest)

AD_imp <- hist_data_imp %>%   # Anderson-Darling test
  group_by(gSSMU, sam.sign) %>%
  summarise(
    n       = n(),
    p_value = ad.test(LnSurvey)$p.value,
    A       = ad.test(LnSurvey)$statistic
  )

# shapiro_results2 <- hist_data_imp %>%  # shapiro.test has a hard limit of 5000 samples
#   group_by(gSSMU, sam.sign) %>%
#   summarise(
#     n = n(),
#     p_value = shapiro.test(LnSurvey)$p.value,
#     W       = shapiro.test(LnSurvey)$statistic
#   )


# log(biomass)
p1b <- ggplot(p1_data, aes(x = log(value), y = after_stat(count / sum(count)))) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = log(1e6), linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_grid(sam_lab ~ gSSMU_lab) +
  coord_cartesian(ylim = c(0, 0.15)) +
  labs(
    title    = "Imputed summer krill biomass (log scale) — Watters et al. (2020) model",
    x        = "ln(Biomass in tonnes)",
    y        = "Relative frequency"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p1b

ggsave(file.path(out.dir, "freq_dist_summer_biomass_log.png"),
       p1b, width = 8, height = 6, dpi = 300)


## 8. Conmparison LogLKB from observed and imputed data
library(MASS)

# Observed data
LKB.obs <- summer_obs %>% dplyr::select(gSSMU, sam.sign, survey) 
LKB.imp <- p1_data %>% dplyr::select(gSSMU_lab, sam_lab, value)
colnames(LKB.obs) <- c("gSSMU", "sam.sign", "biomass")
colnames(LKB.imp) <- c("gSSMU", "sam.sign", "biomass")

# --- Harmonize factor levels before combining ---
LKB.obs <- LKB.obs %>%
  mutate(
    gSSMU    = case_when(
      gSSMU %in% c("1", 1) ~ "gSSMU 1",
      gSSMU %in% c("2", 2) ~ "gSSMU 2"
    ),
    sam.sign = case_when(
      sam.sign %in% c("Neg", "SAM Negative", "negative") ~ "Negative",
      sam.sign %in% c("Pos", "SAM Positive", "positive") ~ "Positive"
    )
  )

LKB.imp <- LKB.imp %>%
  mutate(
    gSSMU    = case_when(
      str_detect(gSSMU, "1") ~ "gSSMU 1",
      str_detect(gSSMU, "2") ~ "gSSMU 2"
    ),
    sam.sign = case_when(
      str_detect(sam.sign, regex("neg", ignore_case = TRUE)) ~ "Negative",
      str_detect(sam.sign, regex("pos", ignore_case = TRUE)) ~ "Positive"
    )
  )


# --- Combine datasets ---
raw_all <- bind_rows(
  LKB.obs %>% mutate(dataset = "LKB.obs"),
  LKB.imp    %>% mutate(dataset = "LKB.imp")
) %>%
  filter(biomass > 0 & !is.na(biomass))

# --- Fit log-normal per dataset × gSSMU × sam.sign ---
fit_lnorm <- function(x) {
  x <- x[x > 0 & !is.na(x)]
  n <- length(x)
  if (n < 3) return(tibble(meanlog = NA_real_, sdlog = NA_real_, n = n))
  fit <- fitdistr(x, "lognormal")
  tibble(meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"], n = n)
}

params_all <- raw_all %>%
  group_by(dataset, gSSMU, sam.sign) %>%
  group_modify(~ fit_lnorm(.x$biomass)) %>%
  ungroup()

# --- Generate fitted curves (normal density on log scale, per facet group) ---
curve_data <- params_all %>%
  filter(!is.na(meanlog),
         gSSMU %in% c("gSSMU 1", "gSSMU 2")) %>%
  rowwise() %>%
  mutate(
    x = list(seq(meanlog - 3.5 * sdlog,
                 meanlog + 3.5 * sdlog,
                 length.out = 300)),
    y = list(dnorm(x, mean = meanlog, sd = sdlog))
  ) %>%
  unnest(c(x, y))

# --- Precompute per-panel per-dataset proportions ---
hist_data <- raw_all %>%
  filter(gSSMU %in% c("gSSMU 1", "gSSMU 2")) %>%
  mutate(bin = floor(log(biomass) / 0.5) * 0.5) %>%
  group_by(gSSMU, sam.sign, dataset, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(gSSMU, sam.sign, dataset) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# --- Precompute means per panel × dataset ---
mean_lines <- raw_all %>%
  filter(gSSMU %in% c("gSSMU 1", "gSSMU 2")) %>%
  group_by(gSSMU, sam.sign, dataset) %>%
  summarise(mean_log = mean(log(biomass), na.rm = TRUE), .groups = "drop")

# --- Plot ---
hist.LKB <- ggplot() +
  geom_col(
    data = hist_data,
    aes(x = bin + 0.25, y = prop, fill = dataset),
    width = 0.5, alpha = 0.5, position = "identity"
  ) +
  geom_line(
    data = curve_data %>% filter(gSSMU %in% c("gSSMU 1", "gSSMU 2")),
    aes(x = x, y = y * 0.5, color = dataset),
    linewidth = 0.9) +
  geom_vline(xintercept = log(1e6),
             linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_vline(
    data = mean_lines,
    aes(xintercept = mean_log, color = dataset),
    linetype = "solid", linewidth = 1) +
  facet_grid(sam.sign ~ gSSMU,
             labeller = labeller(
               gSSMU    = c("gSSMU 1" = "gSSMU 1", "gSSMU 2" = "gSSMU 2"),
               sam.sign = c("Negative" = "SAM Negative", "Positive" = "SAM Positive")
             )) +
  scale_x_continuous(limits = c(x_min, x_max), breaks = x_ll, labels = x_ll) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_fill_manual(values  = c("LKB.obs" = "#b2182b", "LKB.imp" = "#2166ac")) +
  scale_color_manual(values = c("LKB.obs" = "#b2182b", "LKB.imp" = "#2166ac")) +
  labs(
    title = "Krill biomass distribution — observed vs imputed (log scale)",
    x     = "ln(Biomass in tonnes)",
    y     = "Relative frequency",
    fill  = "Dataset", color = "Dataset") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom"
  ); hist.LKB

out.dir <- "./A_Results_Suppl/Model01_filtered/"
ggsave(file.path(out.dir, "hist_LKB_obs_imp.png"),
       hist.LKB,  width=8, height=6, dpi=300)




### We estimate LKB imputed data using Equation (1) and (2) (Waters et al. 2020)

## Equation (2)
# K = U(0.1k, 10k); k = mean(log(LKB)), by SAM sign
# sigma = U(0.1s, 10s); s = sd(log(LKB)), by SAM sign

# First, creating mean(k) and sd(k) matrix
# meanlogsummer <- matrix(
#   nrow = 2,
#   ncol = 2,
#   byrow = TRUE,
#   dimnames = list(
#     c("gSSMU1", "gSSMU2"),
#     c("Neg", "Pos")
#   )
# )
# meanlogsummer[1,] <- Table.01$mu[1:2] 
# meanlogsummer[2,] <- Table.01$mu[3:4] 
# 
# sdlogsummer <- matrix(
#   nrow = 2,
#   ncol = 2,
#   byrow = TRUE,
#   dimnames = list(
#     c("gSSMU1", "gSSMU2"),
#     c("Neg", "Pos")
#   )
# )
# sdlogsummer[1,] <- Table.01$sd[1:2] 
# sdlogsummer[2,] <- Table.01$sd[3:4] 
# 
# # Produce random imputed LKB values for each set o condition (gSSMU x SAM)
# set.seed(123)  # optional, for reproducibility
# 
# n_sim <- 100000
# 
# # meanlogsummer and sdlogsummer are 2x2 matrices:
# # rows: gSSMU1, gSSMU2
# # cols: Neg, Pos
# 
# # log limits
# lower_log <- log(10000)
# upper_log <- log(100000000)
# 
# # prepare empty matrix for simulated LKB
# sim_mat_log <- matrix(NA_real_, nrow = n_sim, ncol = 4)
# colnames(sim_mat_log) <- c("g1.Neg", "g1.Pos", "g2.Neg", "g2.Pos")
# 
# # optional: store mu.ls and sd.ls if you want to keep them
# mu_mat <- sim_mat_log
# sd_mat <- sim_mat_log
# 
# col_index <- 1
# for (i in 1:2) {          # gSSMU: 1,2
#   for (j in 1:2) {        # sam.sign: Neg, Pos
#     # draw mu.ls and sd.ls from uniform distributions
#     mu_ls <- runif(
#       n_sim,
#       min = 0.1 * meanlogsummer[i, j],
#       max = 10  * meanlogsummer[i, j]
#     )
#     sd_ls <- runif(
#       n_sim,
#       min = 0.1 * sdlogsummer[i, j],
#       max = 10  * sdlogsummer[i, j]
#     )
#     
#     # normal draws on log scale
#     log_vals <- rnorm(
#       n    = n_sim,
#       mean = mu_ls,
#       sd   = sd_ls
#     )
#     
#     # truncate to [log(10000), log(1e8)]
#     log_vals <- pmin(pmax(log_vals, lower_log), upper_log)
#     
#     # store
#     sim_mat_log[, col_index] <- log_vals
#     mu_mat[, col_index]      <- mu_ls
#     sd_mat[, col_index]      <- sd_ls
#     
#     col_index <- col_index + 1
#   }
# }
# 
# 
# # back to original (normal) scale
# LKB.imp <- as.data.frame(exp(sim_mat_log))
# head(LKB.imp)
# 
# # Plot histogram of imputed mean-log(LKB):
# 
# # Transform data to be used:
# sim_long <- sim_mat_log %>%
#   as.data.frame() %>%
#   mutate(sim_id = row_number()) %>%
#   pivot_longer(
#     cols = c(g1.Neg, g1.Pos, g2.Neg, g2.Pos),
#     names_to = c("gSSMU", "sam.sign"),
#     names_sep = "\\."
#   ) %>%
#   rename(LnSurvey = value) %>%
#   mutate(
#     gSSMU    = factor(gSSMU, levels = c("g1", "g2"), labels = c(1, 2)),
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
#   )
# 
# #  1. Data for histograms
# hist_imp_data <- sim_long %>%
#   mutate(
#     gSSMU    = factor(gSSMU, levels = c(1, 2)),
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
#   )
# 
# # Test for log-normal distribution
# # AS sample size > 5000, we need to slice the data:
# shapiro_imp <- hist_imp_data %>%
#   group_by(gSSMU, sam.sign) %>%
#   group_modify(~{
#     x       <- .x$LnSurvey
#     n_total <- length(x)
#     
#     # limit to 5000 observations per group
#     if (n_total > 5000) {
#       x <- sample(x, 5000)
#     }
#     
#     st <- shapiro.test(x)
#     
#     tibble(
#       n_total = n_total,
#       n_used  = length(x),
#       W       = unname(st$statistic),
#       p_value = st$p.value
#     )
#   }) %>%
#   ungroup()
# 
# 
# # 2. Data for fitted normal curves, built only from params_ln
# # Imputed LKB 
# Table.02 <- sim_long %>%
#   filter(!is.na(LnSurvey),
#          gSSMU %in% c(1, 2),
#          sam.sign %in% c("Neg", "Pos")) %>%
#   group_by(gSSMU,sam.sign) %>%
#   summarise(
#     n       = n(),
#     mu  = mean(LnSurvey),
#     sd  = sd(LnSurvey),
#     .groups = "drop"
#   )
# tab <- knitr::kable(Table.02, digits = 3, format = "html") |>
#   kableExtra::kable_styling(
#     bootstrap_options = c("striped", "condensed"),
#     full_width = FALSE,
#     font_size = 14
#   )
# 
# # Best: save to PDF (vector, publication-grade)
# kableExtra::save_kable(
#   tab,
#   file = file.path(script_path, "Table_S2.pdf")
# )
# # Or save to PNG (raster)
# kableExtra::save_kable(
#   tab,
#   file = file.path(script_path, "Table_S2.png"),
#   density = 300
# )
# 
# curve_data_imp <- Table.02 %>%
#   mutate(
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
#     gSSMU    = factor(gSSMU)
#   ) %>%
#   rowwise() %>%
#   mutate(
#     x = list(seq(mu - 4*sd, mu + 4*sd, length.out = 200))
#   ) %>%
#   unnest(x) %>%
#   mutate(
#     density = dnorm(x, mean = mu, sd = sd)
#   ) %>%
#   ungroup()
# 
# # 3. Plot: histogram of data + normal line from params_ln
# x_min <- floor(min(hist_imp_data$LnSurvey, na.rm = TRUE))
# x_max <- ceiling(max(hist_imp_data$LnSurvey, na.rm = TRUE))
# x_ll <- seq(x_min, x_max, 1)
# imputed_survey_hist <- ggplot(hist_imp_data, aes(x = LnSurvey)) +
#   geom_histogram(aes(y = after_stat(density)),
#                  binwidth = 0.5,
#                  color = "black",
#                  fill  = "grey80") +
#   geom_line(data = curve_data_imp,
#             aes(x = x, y = density),
#             linewidth = 1) +
#   facet_grid(sam.sign ~ gSSMU) +
#   scale_x_continuous(limits = c(x_min, x_max),
#                      breaks = x_ll, 
#                      labels = x_ll) +
#   labs(x = "LnSurvey", y = "Density") +
#   theme_minimal()
# imputed_survey_hist # Figure S3


### 8. Estimating Frequency distributions for imputed mean-log(LKB) - variation####

### We estimate LKB imputed data using only Equation (1)Waters et al. 2020)
# We use only mean-log(lkb) data; no Uniform distribution

# First, creating mean(k) and sd(k) matrix
# meanlogsummer <- matrix(
#   nrow = 2,
#   ncol = 2,
#   byrow = TRUE,
#   dimnames = list(
#     c("gSSMU1", "gSSMU2"),
#     c("Neg", "Pos")
#   )
# )
# meanlogsummer[1,] <- Table.01$mu[1:2] 
# meanlogsummer[2,] <- Table.01$mu[3:4] 
# 
# sdlogsummer <- matrix(
#   nrow = 2,
#   ncol = 2,
#   byrow = TRUE,
#   dimnames = list(
#     c("gSSMU1", "gSSMU2"),
#     c("Neg", "Pos")
#   )
# )
# sdlogsummer[1,] <- Table.01$sd[1:2] 
# sdlogsummer[2,] <- Table.01$sd[3:4] 
# 
# # Produce random imputed LKB values for each set o condition (gSSMU x SAM)
# set.seed(123)  # optional, for reproducibility
# 
# n_sim <- 1000
# 
# # meanlogsummer and sdlogsummer are 2x2 matrices:
# # rows: gSSMU1, gSSMU2
# # cols: Neg, Pos
# 
# # log limits
# lower_log <- log(10000)
# upper_log <- log(100000000)
# 
# # prepare empty matrix for simulated LKB
# sim_mat_log2 <- matrix(NA_real_, nrow = n_sim, ncol = 4)
# colnames(sim_mat_log2) <- c("g1.Neg", "g1.Pos", "g2.Neg", "g2.Pos")
# 
# # optional: store mu.ls and sd.ls if you want to keep them
# mu_mat <- sim_mat_log
# sd_mat <- sim_mat_log
# 
# col_index <- 1
# for (i in 1:2) {          # gSSMU: 1,2
#   for (j in 1:2) {        # sam.sign: Neg, Pos
#     # draw mu.ls and sd.ls from meanlogsummer / sdlogsummer
#     mu_ls <- meanlogsummer[i,j]
#     sd_ls <- sdlogsummer[i,j]
#     
#     # normal draws on log scale
#     log_vals <- rnorm(
#       n    = n_sim,
#       mean = mu_ls,
#       sd   = sd_ls
#     )
#     
#     # truncate to [log(10000), log(1e8)]
#     log_vals <- pmin(pmax(log_vals, lower_log), upper_log)
#     
#     # store
#     sim_mat_log2[, col_index] <- log_vals
#     mu_mat[, col_index]      <- mu_ls
#     sd_mat[, col_index]      <- sd_ls
#     
#     col_index <- col_index + 1
#   }
# }
# 
# 
# # back to original (normal) scale
# LKB.imp.mean <- as.data.frame(exp(sim_mat_log2))
# head(LKB.imp.mean)
# 
# # Plot histogram of imputed mean-log(LKB):
# 
# # Transform data to be used:
# sim_long2 <- sim_mat_log2 %>%
#   as.data.frame() %>%
#   mutate(sim_id = row_number()) %>%
#   pivot_longer(
#     cols = c(g1.Neg, g1.Pos, g2.Neg, g2.Pos),
#     names_to = c("gSSMU", "sam.sign"),
#     names_sep = "\\."
#   ) %>%
#   rename(LnSurvey = value) %>%
#   mutate(
#     gSSMU    = factor(gSSMU, levels = c("g1", "g2"), labels = c(1, 2)),
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
#   )
# 
# #  1. Data for histograms
# hist_imp_data2 <- sim_long2 %>%
#   mutate(
#     gSSMU    = factor(gSSMU, levels = c(1, 2)),
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos"))
#   )
# 
# # 2. Data for fitted normal curves, built only from params_ln
# # Imputed LKB 
# Table.03 <- sim_long2 %>%
#   filter(!is.na(LnSurvey),
#          gSSMU %in% c(1, 2),
#          sam.sign %in% c("Neg", "Pos")) %>%
#   group_by(gSSMU,sam.sign) %>%
#   summarise(
#     n       = n(),
#     mu  = mean(LnSurvey),
#     sd  = sd(LnSurvey),
#     .groups = "drop"
#   )
# knitr::kable(Table.03, digits = 3)
# 
# curve_data_imp2 <- Table.03 %>%
#   mutate(
#     sam.sign = factor(sam.sign, levels = c("Neg", "Pos")),
#     gSSMU    = factor(gSSMU)
#   ) %>%
#   rowwise() %>%
#   mutate(
#     x = list(seq(mu - 4*sd, mu + 4*sd, length.out = 200))
#   ) %>%
#   unnest(x) %>%
#   mutate(
#     density = dnorm(x, mean = mu, sd = sd)
#   ) %>%
#   ungroup()
# 
# # 3. Plot: histogram of data + normal line from params_ln
# x_min <- floor(min(hist_imp_data2$LnSurvey, na.rm = TRUE))
# x_max <- ceiling(max(hist_imp_data2$LnSurvey, na.rm = TRUE))
# x_ll <- seq(x_min, x_max, 1)
# imputed_survey_hist.mean <- ggplot(hist_imp_data2, aes(x = LnSurvey)) +
#   geom_histogram(aes(y = after_stat(density)),
#                  binwidth = 0.5,
#                  color = "black",
#                  fill  = "grey80") +
#   geom_line(data = curve_data_imp2,
#             aes(x = x, y = density),
#             linewidth = 1) +
#   facet_grid(sam.sign ~ gSSMU) +
#   scale_x_continuous(limits = c(x_min, x_max),
#                      breaks = x_ll, 
#                      labels = x_ll) +
#   labs(x = "LnSurvey", y = "Density") +
#   theme_minimal()
# imputed_survey_hist.mean
# 





### 9. Predator's krill consumption ####

Table.05 <- read.csv("./A_Results_Suppl/Baleen Whale Abundance.csv",
                     header=TRUE, fileEncoding = "Windows-1252")


Predators_summer <- read.csv("./A_Results_Suppl/Summer_Krill_req.csv",header=TRUE)
Predators_winter <- read.csv("./A_Results_Suppl/Winter_Krill_req.csv",header=TRUE)

Table.06 <- Predators_summer %>% select(-Fish_Qsummer_g, -Whale_Qsummer_g, -Penguin_Qsummer_g, -Reference)
# knitr::kable(Table.03, digits = 0)

Table.07 <- Predators_winter %>% select(-Fish_Qwinter_g, -Whale_Qwinter_g, -Penguin_Qwinter_g, -Reference)
# knitr::kable(Table.04, digits = 0)


