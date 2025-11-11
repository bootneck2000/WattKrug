# --- Clean & aggregate survey data ---
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
survey <- survey_low_nmi %>%
  select(Year, Leg, gSSMU, biomass, nmi.count)
survey$season<-ifelse(survey$Year<2012,"S","W")
survey <- survey %>% rename(cal.yr = Year)  # to match rest of script

survey <- survey %>%
  mutate(
    gSSMU  = as.character(gSSMU),
    season = toupper(season)
  ) %>%
  group_by(cal.yr, gSSMU, season) %>%
  summarise(survey = mean(biomass, na.rm = TRUE), .groups = "drop")
rm(percentiles, survey_filtered, survey_low_nmi)
