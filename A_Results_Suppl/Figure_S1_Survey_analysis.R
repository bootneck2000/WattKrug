##### KRILL SURVEY AND FISHERY DATA

library(dplyr)
library(ggplot2)
library(patchwork)

## Analysis of acoustic coverage per gSSSMU, leg and year

# krill survey biomass
survey<-read.csv("./Supplementary Files/krillsurveywithJoinville.csv",header=TRUE,stringsAsFactors = FALSE)


# --- Prepare data (gSSMU = 2) -----------------------------------------------
g <- 2
by_year_leg_2 <- survey %>%
  filter(gSSMU == g) %>%
  group_by(Year, Leg) %>%
  summarise(nmi_total = sum(nmi.count, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(nmi_total))

# Observed years in this gSSMU (to keep X relevant)
years_2 <- by_year_leg_2 %>% pull(Year) %>% unique() %>% sort()

# 10th percentile threshold across Year×Leg (within this gSSMU)
q10_2 <- quantile(by_year_leg_2$nmi_total, 0.10, na.rm = TRUE)

# Histogram x breaks (binwidth = 10) tailored to this gSSMU
max_2 <- max(by_year_leg_2$nmi_total, na.rm = TRUE)
max_2_round <- round(max_2,0) # ceiling(max_2 / 10) * 10
hist_breaks_2 <- seq(0, max_2_round, by = 50)

# --- (A) Bar plot: Year × Leg -----------------------------------------------
p1_2 <- ggplot(by_year_leg_2,
               aes(x = factor(Year, levels = years_2),
                   y = nmi_total,
                   fill = factor(Leg))) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = round(q10_2,0), color = "black",  
             linewidth = 0.5) + #, linetype = "dashed") +
  labs(
    title = "gSSMU 2 — Survey effort",
    x = "Year", y = "Total transect length (nm)", fill = "Leg"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- (B) Coverage histogram: Year×Leg totals, binwidth = 10 -----------------
p2_2 <- ggplot(by_year_leg_2, aes(x = nmi_total)) +
  geom_histogram(binwidth = 10, boundary = 0, closed = "left") +
  geom_vline(xintercept = round(q10_2,0), linetype = "dashed") +
  scale_x_continuous(breaks = hist_breaks_2, limits = c(0, max_2_round)) +
  labs(
    title = "Total transect length",
    x = "Sum of transect length (nm)", y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_g2 <- p1_2 / p2_2 + plot_layout(heights = c(2, 1))
# Show
fig_g2


# --- Prepare data (gSSMU = 1) -----------------------------------------------
g <- 1
by_year_leg_1 <- survey %>%
  filter(gSSMU == g) %>%
  group_by(Year, Leg) %>%
  summarise(nmi_total = sum(nmi.count, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(nmi_total))

# Observed years in this gSSMU (to keep X relevant)
years_1 <- by_year_leg_1 %>% pull(Year) %>% unique() %>% sort()

# 10th percentile threshold across Year×Leg (within this gSSMU)
q10_1 <- quantile(by_year_leg_1$nmi_total, 0.10, na.rm = TRUE)

# Histogram x breaks (binwidth = 10) tailored to this gSSMU
max_1 <- max(by_year_leg_1$nmi_total, na.rm = TRUE)
max_1_round <- round(max_1,0) # ceiling(max_1 / 10) * 10
hist_breaks_1 <- seq(0, max_1_round, by = 50)

# --- (A) Bar plot: Year × Leg -----------------------------------------------
p1_1 <- ggplot(by_year_leg_1,
               aes(x = factor(Year, levels = years_1),
                   y = nmi_total,
                   fill = factor(Leg))) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = round(q10_1,0), color = "black",  
             linewidth = 0.5) + #, linetype = "dashed") +
  labs(
    title = "gSSMU 1 — Survey effort",
    x = "Year", y = "Total transect length (nm)", fill = "Leg"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- (B) Coverage histogram: Year×Leg totals, binwidth = 10 -----------------
p2_1 <- ggplot(by_year_leg_1, aes(x = nmi_total)) +
  geom_histogram(binwidth = 10, boundary = 0, closed = "left") +
  geom_vline(xintercept = round(q10_1,0), linetype = "dashed") +
  scale_x_continuous(breaks = hist_breaks_1, limits = c(0, max_1_round)) +
  labs(
    title = "Total transect length",
    x = "Sum of transect length (nm)", y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_g1 <- p1_1 / p2_1 + plot_layout(heights = c(2, 1))
# Show
fig_g1




