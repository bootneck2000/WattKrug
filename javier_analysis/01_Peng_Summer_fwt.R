##### Review Penguin indices
##### generate SUMMER INDICES
##### fledge weight (fwt) 

## Watters Code ## Do not change
# bigger indicates better summer

fwt <- read.csv("./Supplementary Files/fweight.csv", header=TRUE, stringsAsFactors = FALSE) # Fledgling mass
fwt<-tapply(fwt$WT,list(fwt$YEAR,fwt$PROJECT,fwt$SPECIES),mean)
fwt<-data.frame(YEAR=rep(dimnames(fwt)[[1]],dim(fwt)[2]*dim(fwt)[3]),
                PROJECT=rep(rep(dimnames(fwt)[[2]],each=dim(fwt)[1]),dim(fwt)[3]),
                SPECIES=rep(dimnames(fwt)[[3]],each=dim(fwt)[1]*dim(fwt)[2]),
                fwt=c(fwt),stringsAsFactors = FALSE)
### They estimated the mean FWT per year. And from that they estimated the mean FWT for all years (the mean of the mean).
fwt$matchme<-paste(fwt$PROJECT,fwt$SPECIES,sep="|")
tt<-tapply(fwt$fwt,list(fwt$matchme),mean,na.rm=TRUE)
ttt<-tapply(fwt$fwt,list(fwt$matchme),sd,na.rm=TRUE)
mean.fwt<-tt[match(fwt$matchme,names(tt))]
sd.fwt<-ttt[match(fwt$matchme,names(ttt))]
fwt$std.mean.fwt<-(fwt$fwt-mean.fwt)/sd.fwt  # INDEX
fwt<-fwt[,-c(4:5)]
#omits<-(fwt$SPECIES=="ADPE"&fwt$PROJECT=="CS")|(fwt$SPECIES=="CHPE"&fwt$PROJECT=="COPA")
#fwt<-fwt[!omits,]
names(fwt)[4]<-"index"
fwt$param=rep("FWT",dim(fwt)[1])
fwt$season=rep("S",dim(fwt)[1])
# make stuff reference the correct "calendar year" for matching up with krill survey and catch data
# summer indices are relevant to the second year in the split-season designation
fwt$cal.yr<-as.numeric(substr(fwt$YEAR,1,4))+1
#print(str(fwt))

## End of Watters Code ##


# ==================== ANALYSIS  ====================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(officer)
  library(stringr)
  library(rlang)
})
suppressWarnings(RNGkind("default"))  # silence harmless RNG warnings on some R builds

# -------------------- Preconditions & setup ------------------------------
stopifnot(exists("fwt"), is.data.frame(fwt))

# --- Output locations ---------------------------
out_dir <- "./javier_analysis/review_FWT"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Data Prep ----------------------------------------
# Read the same source path used by the Watters code
ade <- read.csv("./Supplementary Files/fweight.csv", header=TRUE, stringsAsFactors = FALSE)

ade <- ade |>
  mutate(
    DATE        = suppressWarnings(lubridate::mdy(DATE)),
    cal.yr      = as.integer(substr(YEAR, 1, 4))+1,
    julian_days = lubridate::yday(DATE)
  )

# ade_m <- ade %>% select(cal.yr, julian_days, PROJECT, SPECIES, WT, AGE)

# Step 1: YEAR x PROJECT x SPECIES means
step1 <- ade %>%
  filter(!is.na(WT)) %>%
  group_by(cal.yr, PROJECT, SPECIES) %>%
  summarise(
    WTmean   = mean(WT, na.rm = TRUE),
    WTsd     = sd(WT, na.rm = TRUE),
    AGEmean   = mean(AGE, na.rm = TRUE),
    AGEsd     = sd(AGE, na.rm = TRUE),
    DAYmean  = mean(julian_days, na.rm = TRUE),
    DAYsd = sd(julian_days, na.rm = TRUE), 
    .groups  = "drop"
  )

# Step 2: Within PROJECT x SPECIES: global mean/sd and day mean
step2 <- step1 %>%
  group_by(PROJECT, SPECIES) %>%
  summarise(
    mean_WT = mean(WTmean, na.rm = TRUE),
    sd_WT   = sd(WTmean, na.rm = TRUE),
    mean_AGE = mean(AGEmean, na.rm = TRUE),
    sd_AGE   = sd(AGEmean, na.rm = TRUE),
    mean_day= mean(DAYmean, na.rm = TRUE),
    sd_day   = sd(DAYmean, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Join & indices (guard against sd_WT == 0)
ADE_sum <- step1 %>%
  left_join(step2, by = c("PROJECT","SPECIES")) %>%
  mutate(
    index     = ifelse(sd_WT > 0, (WTmean - mean_WT)/sd_WT, NA_real_),
    AGE_index = ifelse(sd_AGE > 0, (AGEmean - mean_AGE)/sd_AGE, NA_real_),
    DAY_index = ifelse(sd_day > 0, (DAYmean - mean_day)/sd_day, NA_real_)
  ) %>%
  arrange(SPECIES, PROJECT, cal.yr)

# -------------------- CORRELATIONS ------------------------------
# Robust correlation-by-group helper with status flag
source(file.path("./javier_analysis/corr_tbl_ts.R"))

# choose candidate variables (means)
cand <- c("cal.yr", "index", "AGE_index", "DAY_index",
                    "WTmean","AGEmean","DAYmean")
cand <- intersect(cand, names(ADE_sum))
pairs <- combn(cand, 2, simplify = FALSE)

run_corr_pair <- function(df, nm_x, nm_y) {
  x_sym <- sym(nm_x); y_sym <- sym(nm_y)
  corr_tbl_ts(df, !!x_sym, !!y_sym) %>% mutate(x_var=nm_x, y_var=nm_y, .before=1)
}

all_corrs   <- map_dfr(pairs, ~ run_corr_pair(ADE_sum,   .x[1], .x[2]))

# Save
write_csv(all_corrs, file.path(out_dir, "FWT_all_combinations_corr.csv"))

# -------------------- PLOTTING SIGNIFICANT PAIRS ----------
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "-", x)

# Make one plot per (x_var, y_var) pair that has ANY significant group
# Include ALL SPECIES x PROJECT in the plot; mark significant panels with a ★
save_pair_plots_all_groups <- function(corr_tbl, data_df) {
  
  # which pairs have at least one significant result?
  sig_pairs <- corr_tbl %>%
    filter(lm_significant == "significant") %>%
    distinct(x_var, y_var)
  
  if (nrow(sig_pairs) == 0L) {
    message("No significant pairs found in FWT")
    return(invisible(NULL))
  }
  
  pwalk(sig_pairs, function(x_var, y_var) {
    
    # build a panel-level significance flag for this pair
    sig_flags <- corr_tbl %>%
      filter(x_var == !!x_var, y_var == !!y_var, lm_R2 < 1) %>%
      mutate(sig_flag = lm_significant == "significant") %>%
      select(SPECIES, PROJECT, sig_flag) %>%
      distinct()
    
    # Filter data to rows that have both variables available
    dat <- data_df %>%
      filter(is.finite(.data[[x_var]]), is.finite(.data[[y_var]])) %>%
      # keep only groups that are present in corr table (optional, but tight)
      inner_join(sig_flags %>% select(SPECIES, PROJECT) %>% distinct(),
                 by = c("SPECIES", "PROJECT")) %>%
      left_join(sig_flags, by = c("SPECIES", "PROJECT"))
    
    if (nrow(dat) == 0L) {
      message(sprintf("Skipping — %s vs %s: no data rows.", y_var, x_var))
      return(invisible(NULL))
    }
    
    # A tiny df to add facet-specific star in top-right corner
    # We'll place label at the 97.5th percentile of each panel to avoid clipping.
    label_df <- dat %>%
      group_by(SPECIES, PROJECT) %>%
      summarise(
        x_pos = suppressWarnings(stats::quantile(.data[[x_var]], 0.97, na.rm = TRUE)),
        y_pos = suppressWarnings(stats::quantile(.data[[y_var]], 0.97, na.rm = TRUE)),
        sig_flag = any(sig_flag %in% TRUE),
        .groups = "drop"
      ) %>%
      mutate(label = ifelse(sig_flag, "\u2605", "")) # ★
    
    p <- ggplot(dat, aes(x = .data[[x_var]], y = .data[[y_var]], color = SPECIES)) +
      geom_point(alpha = 0.9, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.25) +
      geom_smooth(method = "loess", se = TRUE, span = 0.75, color = "lightgreen", alpha = 0.15) +
      geom_text(
        data = label_df,
        aes(x = x_pos, y = y_pos, label = label),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 5,
        show.legend = FALSE
      ) +
      facet_grid(SPECIES ~ PROJECT, scales = "free_y") +
      labs(
        title = sprintf("FWT — %s vs %s (all groups; ★ = significant)", y_var, x_var),
        x = x_var, y = y_var
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold")
      )
    
    fpath <- file.path(out_dir, paste0("FWT_", sanitize(y_var), "_", sanitize(x_var), ".png"))
    ggsave(fpath, p, width = 12, height = 8, dpi = 300, bg = "white")
    message("Saved: ", fpath)
  })
}

save_pair_plots_all_groups(all_corrs, ADE_sum)


## --- Although NO significant, I added the cal.yr versus indexes
save_pair_plots_all_groups <- function(corr_tbl, data_df) {
  
  # which pairs have at least one significant result?
  sig_pairs <- corr_tbl %>%
    filter(x_var == "cal.yr") %>%
    distinct(x_var, y_var)
  
  if (nrow(sig_pairs) == 0L) {
    message("No pairs found in FWT")
    return(invisible(NULL))
  }
  
  pwalk(sig_pairs, function(x_var, y_var) {
    
    # build a panel-level significance flag for this pair
    sig_flags <- corr_tbl %>%
      filter(x_var == !!x_var, y_var == !!y_var, lm_R2 < 1) %>%
      mutate(sig_flag = lm_significant == "significant") %>%
      select(SPECIES, PROJECT, sig_flag) %>%
      distinct()
    
    # Filter data to rows that have both variables available
    dat <- data_df %>%
      filter(is.finite(.data[[x_var]]), is.finite(.data[[y_var]])) %>%
      # keep only groups that are present in corr table (optional, but tight)
      inner_join(sig_flags %>% select(SPECIES, PROJECT) %>% distinct(),
                 by = c("SPECIES", "PROJECT")) %>%
      left_join(sig_flags, by = c("SPECIES", "PROJECT"))
    
    if (nrow(dat) == 0L) {
      message(sprintf("Skipping — %s vs %s: no data rows.", y_var, x_var))
      return(invisible(NULL))
    }
    
    # A tiny df to add facet-specific star in top-right corner
    # We'll place label at the 97.5th percentile of each panel to avoid clipping.
    label_df <- dat %>%
      group_by(SPECIES, PROJECT) %>%
      summarise(
        x_pos = suppressWarnings(stats::quantile(.data[[x_var]], 0.97, na.rm = TRUE)),
        y_pos = suppressWarnings(stats::quantile(.data[[y_var]], 0.97, na.rm = TRUE)),
        sig_flag = any(sig_flag %in% TRUE),
        .groups = "drop"
      ) %>%
      mutate(label = ifelse(sig_flag, "\u2605", "")) # ★
    
    p <- ggplot(dat, aes(x = .data[[x_var]], y = .data[[y_var]], color = SPECIES)) +
      geom_point(alpha = 0.9, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.25) +
      geom_smooth(method = "loess", se = TRUE, span = 0.75, color = "lightgreen", alpha = 0.15) +
      geom_text(
        data = label_df,
        aes(x = x_pos, y = y_pos, label = label),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 5,
        show.legend = FALSE
      ) +
      facet_grid(SPECIES ~ PROJECT, scales = "free_y") +
      labs(
        title = sprintf("FWT — %s vs %s (all groups; ★ = significant)", y_var, x_var),
        x = x_var, y = y_var
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold")
      )
    
    fpath <- file.path(out_dir, paste0("FWT_", sanitize(y_var), "_", sanitize(x_var), ".png"))
    ggsave(fpath, p, width = 12, height = 8, dpi = 300, bg = "white")
    message("Saved: ", fpath)
  })
}

save_pair_plots_all_groups(all_corrs, ADE_sum)


# -------------------- Effect of Years on Y-value --------------------
# For each significant (x_var, y_var) within each SPECIES × PROJECT,
# compare mean(y_var) for first and seccond half of the study period
compute_half_means <- function(corr_df, data_df) {
  sig_pairs <- corr_df %>%
    filter(lm_significant == "significant", lm_R2 < 1) %>%
    select(SPECIES, PROJECT, x_var, y_var, n_ok, lm_slope, lm_p) %>%
    distinct()
  
  if (!nrow(sig_pairs)) return(tibble())
  
  pmap_dfr(sig_pairs, function(SPECIES, PROJECT, x_var, y_var, n_ok, lm_slope, lm_p) {
    dat <- data_df %>%
      filter(SPECIES == !!SPECIES, PROJECT == !!PROJECT) %>%
      filter(is.finite(.data[[x_var]]), is.finite(.data[[y_var]]))
    
    if (!nrow(dat)) return(tibble())
    
    first_year <- min(dat$cal.yr, na.rm = TRUE)
    last_year  <- max(dat$cal.yr, na.rm = TRUE)
    mid_year   <- floor((first_year + last_year)/2)
    
    dat %>%
      mutate(period = ifelse(cal.yr <= mid_year, "first_half", "second_half")) %>%
      group_by(period) %>%
      summarise(
        SPECIES    = first(SPECIES),
        PROJECT    = first(PROJECT),
        x_var      = first(x_var),
        y_var      = first(y_var),
        n_used     = n(),
        mean_y     = mean(.data[[y_var]], na.rm = TRUE),
        sd_y       = sd(.data[[y_var]], na.rm = TRUE),
        first_year = first_year,
        last_year  = last_year,
        mid_year   = mid_year,
        lm_slope   = lm_slope,
        lm_p       = lm_p,
        .groups    = "drop"
      )
  })
}

halves <- compute_half_means(all_corrs, ADE_sum)

# Neat table
halves_neat <- halves %>%
  mutate(
    lm_slope = round(lm_slope, 4),
    lm_p     = signif(lm_p, 3),
    mean_y   = round(mean_y, 3),
    sd_y     = round(sd_y, 3)
  ) %>%
  select(
    SPECIES, PROJECT, x_var, y_var,
    first_year, mid_year, last_year, lm_slope, lm_p,
    period, n_used, mean_y, sd_y
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = c(n_used, mean_y, sd_y),
    names_sep = "_"
  ) %>%
  # rename for readability
  rename(
    n_first     = n_used_first_half,
    n_second    = n_used_second_half,
    mean_y_first  = mean_y_first_half,
    mean_y_second = mean_y_second_half,
    sd_y_first    = sd_y_first_half,
    sd_y_second   = sd_y_second_half
  ) %>%
  arrange(SPECIES, PROJECT, x_var, y_var) %>%
  mutate(Diff_range = mean_y_first - mean_y_second)

# Calculate Relative Difference (Normalized by Mean) and Coefficient of Variation–like Metric
# Summarise by group
DatMean <- ADE_sum %>%
  group_by(SPECIES, PROJECT, mean_WT, sd_WT, mean_AGE, sd_AGE, mean_day, sd_day) %>%
  summarise(
    Index_mean     = mean(index, na.rm = TRUE),
    Index_sd     = sd(index, na.rm = TRUE),
    AGE_Index_mean     = mean(AGE_index, na.rm = TRUE),
    AGE_Index_sd     = sd(AGE_index, na.rm = TRUE),
    DAY_Index_mean = mean(DAY_index, na.rm = TRUE),
    DAY_Index_sd = sd(DAY_index, na.rm = TRUE),
    .groups = "drop"
  )

# Left-join group averages into sig_pairs_mean_y_neat
sig_pairs_with_group_means <- halves_neat %>%
  left_join(DatMean,
            by = c("SPECIES", "PROJECT"))

## Estimate Relative Difference (Normalized by Mean)
#Simple, dimensionless.
#Interpreted as “percentage of the mean.”
#Sensitive if mean ≈ 0.
unique(sig_pairs_with_group_means$y_var)
sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  mutate(
    RD = case_when(
      y_var == "DAY_index" ~ (Diff_range / DAY_Index_mean),
      y_var == "DAYmean"   ~ (Diff_range / mean_day),
      y_var == "WTmean"    ~ (Diff_range / mean_WT),
      y_var == "AGEmean"    ~ (Diff_range / mean_AGE),
      y_var == "AGE_index"    ~ (Diff_range / AGE_Index_mean),
      y_var == "index"     ~ (Diff_range / Index_mean),
      TRUE ~ NA_real_   # fallback for safety
    )
  )

## Estimate Coefficient of Variation–like Metric
#Normalizes by standard deviation of all values.
#Interpreted as “how many SDs apart are these two values.”
sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  mutate(
    CVdiff = case_when(
      y_var == "DAY_index" ~ abs(Diff_range) / DAY_Index_sd,
      y_var == "DAYmean"   ~ abs(Diff_range) / sd_day,
      y_var == "WTmean"    ~ abs(Diff_range) / sd_WT,
      y_var == "AGEmean"    ~ abs(Diff_range) / sd_AGE,
      y_var == "AGE_index"    ~ abs(Diff_range) / AGE_Index_sd,
      y_var == "index"     ~ abs(Diff_range) / Index_sd,
      TRUE ~ NA_real_   # fallback for safety
    )
  )

sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  mutate(RD = round(RD,5)) %>%
  mutate(CVdiff = round(CVdiff,3))

sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  select(-mean_WT, -sd_WT, -mean_day, -sd_day, -Index_mean, -Index_sd, 
         -DAY_Index_mean, -DAY_Index_sd, -RD, -AGE_Index_sd, -AGE_Index_mean,
         -mean_AGE, -sd_AGE)

# Save and pretty-print
write_csv(sig_pairs_with_group_means, file.path(out_dir, "FWT_significant_pairs_mean_y_by_halves.csv"))

# ==================== WORD REPORT (OFFICER) =======================================
docx_path <- file.path(out_dir, "FWT_Correlations_Report.docx")

round_numeric_cols <- function(df, digits = 3) {
  dplyr::mutate(df, dplyr::across(where(is.numeric), ~round(.x, digits)))
}

# Build clean, compact tables for Word
Comparisons <- sig_pairs_with_group_means %>%
  dplyr::select(SPECIES, PROJECT, x_var, y_var, lm_slope, lm_p, Diff_range, CVdiff) %>%
  round_numeric_cols(3)

# Create the Word document
ps_landscape <- prop_section(
  page_size   = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.6, right = 0.6, bottom = 0.6, left = 0.6)
)

doc <- read_docx(path = system.file("template/template.docx", package = "officer"))
doc <- body_set_default_section(doc, ps_landscape)

doc <- body_add_par(doc, "FWT — Correlations & Half-Period Comparison", style = "heading 1")
doc <- body_add_par(doc, "Significant pairs & first vs second period comparison", style = "heading 2")
doc <- body_add_table(doc, value = Comparisons, style = "Normal Table")

# Add the “all groups” pair-plot files if you generated them
doc <- body_add_par(doc, "Significant pairs — facet plots (★ marks significant panels)", style = "heading 2")
for (img in list.files(out_dir, pattern = "^FWT_.*\\.png$", full.names = TRUE)) {
  doc <- body_add_img(doc, src = img, width = 7, height = 4.7)
}

print(doc, target = docx_path)
message("Word report saved at: ", normalizePath(docx_path, winslash = "/", mustWork = FALSE))

