##### Review Penguin indices
##### generate SUMMER INDICES
##### foraging trip duration (td)


## Watters Code:
# smaller indicates better summer (thus need to switch direction of index)
td <- read.csv("./Supplementary Files/tripduration.csv") # Forging trip duration
td<-td[,c(1:3,8)]
# next line is to make trip duration point in same direction as fwt and phs (max td is 59.95 for all trips)
# call this "revtd" for "reversed" trip duration
td[,4]<-60-td[,4]     # reverse
names(td)[4]<-"revtd"
td<-tapply(td$revtd,list(td$YEAR,td$PROJECT,td$SPECIES),mean)
td<-data.frame(YEAR=rep(dimnames(td)[[1]],dim(td)[2]*dim(td)[3]),
               PROJECT=rep(rep(dimnames(td)[[2]],each=dim(td)[1]),dim(td)[3]),
               SPECIES=rep(dimnames(td)[[3]],each=dim(td)[1]*dim(td)[2]),
               revtd=c(td),stringsAsFactors = FALSE)
td$matchme<-paste(td$PROJECT,td$SPECIES,sep="|")
tt<-tapply(td$revtd,list(td$matchme),mean,na.rm=TRUE)
ttt<-tapply(td$revtd,list(td$matchme),sd,na.rm=TRUE)
mean.revtd<-tt[match(td$matchme,names(tt))]
sd.revtd<-ttt[match(td$matchme,names(ttt))]
td$std.revtd<-(td$revtd-mean.revtd)/sd.revtd  # INDEX
td<-td[,-c(4:5)]
names(td)[4]<-"index"
#omits<-(td$SPECIES=="ADPE"&td$PROJECT=="CS")|(td$SPECIES=="CHPE"&td$PROJECT=="COPA")
#td<-td[!omits,]
td$param=rep("REVTD",dim(td)[1])
td$season=rep("S",dim(td)[1])
# summer indices are relevant to the second year in the split-season designation
td$cal.yr<-as.numeric(substr(td$YEAR,1,4))+1

## End of Watters Code ##
plot(td$cal.yr, td$index)

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
stopifnot(exists("td"), is.data.frame(td))

# --- Output locations ---------------------------
out_dir <- "./javier_analysis/review_TD"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Data Prep ----------------------------------------
# Read the same source path used by the Watters code
ade <- read.csv("./Supplementary Files/tripduration.csv") 

ade <- ade |>
  mutate(
    DATE        = suppressWarnings(lubridate::mdy(DATE)),
    cal.yr      = as.integer(substr(YEAR, 1, 4))+1,
    julian_days = lubridate::yday(DATE),
    revtd       = 60 - TRIP_DURATION
  )

# Step 1: YEAR x PROJECT x SPECIES means
step1 <- ade %>%
  filter(!is.na(TRIP_DURATION)) %>%
  group_by(cal.yr, PROJECT, SPECIES) %>%
  summarise(
    DAYseason_mean   = mean(DAY_OF_SEASON, na.rm = TRUE),
    WEEKseason_mean   = mean(WEEK_OF_SEASON, na.rm = TRUE),
    TD_mean  = mean(TRIP_DURATION, na.rm = TRUE),
    REVtd_mean = mean(revtd, na.rm = TRUE),
    .groups  = "drop"
  )

# Step 2: Within PROJECT x SPECIES: global mean/sd and day mean
step2 <- step1 %>%
  group_by(PROJECT, SPECIES) %>%
  summarise(
    mean_DAYs = mean(DAYseason_mean, na.rm = TRUE),
    sd_DAYs   = sd(DAYseason_mean, na.rm = TRUE),
    mean_WKs = mean(WEEKseason_mean, na.rm = TRUE),
    sd_WKs   = sd(WEEKseason_mean, na.rm = TRUE),
    mean_TD_ind = mean(TD_mean, na.rm = TRUE),
    sd_TD_ind   = sd(TD_mean, na.rm = TRUE),
    mean_REVtd = mean(REVtd_mean, na.rm = TRUE), 
    sd_REVtd = sd(REVtd_mean, na.rm = TRUE), 
    .groups = "drop"
  )

# Step 3: Join & indices (guard against sd_WT == 0)
ADE_sum <- step1 %>%
  left_join(step2, by = c("PROJECT","SPECIES")) %>%
  mutate(
    DAYs_index     = ifelse(sd_DAYs > 0, (DAYseason_mean - mean_DAYs)/sd_DAYs, NA_real_),
    WKs_index     = ifelse(sd_WKs > 0, (WEEKseason_mean - mean_WKs)/sd_WKs, NA_real_),
    REVtd_index     = ifelse(sd_REVtd > 0, (REVtd_mean - mean_REVtd)/sd_REVtd, NA_real_),
  ) %>%
  arrange(SPECIES, PROJECT, cal.yr)

# -------------------- CORRELATIONS ------------------------------
# Robust correlation-by-group helper with status flag
source(file.path("./javier_analysis/corr_tbl_ts.R"))

# choose candidate variables (means)
colnames(ADE_sum)
cand <- c("cal.yr", "DAYs_index", "WKs_index",  "REVtd_index", "mean_DAYs", 
          "mean_WKs", "mean_TD_ind", "mean_REVtd")

cand <- intersect(cand, names(ADE_sum))
pairs <- combn(cand, 2, simplify = FALSE)

run_corr_pair <- function(df, nm_x, nm_y) {
  x_sym <- sym(nm_x); y_sym <- sym(nm_y)
  corr_tbl_ts(df, !!x_sym, !!y_sym) %>% mutate(x_var=nm_x, y_var=nm_y, .before=1)
}

all_corrs   <- map_dfr(pairs, ~ run_corr_pair(ADE_sum,   .x[1], .x[2]))

# Save
write_csv(all_corrs, file.path(out_dir, "TD_all_combinations_corr.csv"))

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
    message("No significant pairs found in TD")
    return(invisible(NULL))
  }
  
  pwalk(sig_pairs, function(x_var, y_var) {
    
    # build a panel-level significance flag for this pair
    sig_flags <- corr_tbl %>%
      filter(x_var == !!x_var, y_var == !!y_var, lm_R2 < 0.90) %>%
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
        title = sprintf("PHS — %s vs %s (all groups; ★ = significant)", y_var, x_var),
        x = x_var, y = y_var
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "lightgrey", color = "black"),
        strip.text = element_text(face = "bold")
      )
    
    fpath <- file.path(out_dir, paste0("TD_", sanitize(y_var), "_", sanitize(x_var), ".png"))
    ggsave(fpath, p, width = 12, height = 8, dpi = 300, bg = "white")
    message("Saved: ", fpath)
  })
}

save_pair_plots_all_groups(all_corrs, ADE_sum)

# -------------------- Effect of Years on Y-value --------------------
# For each significant (x_var, y_var) within each SPECIES × PROJECT,
# compare mean(y_var) for first and second half of the study period
compute_half_means <- function(corr_df, data_df) {
  sig_pairs <- corr_df %>%
    filter(lm_significant == "significant", lm_R2 < 0.9) %>%
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
colnames(ADE_sum)
DatMean <- ADE_sum %>%
  group_by(SPECIES, PROJECT, mean_DAYs, sd_DAYs, mean_WKs,  sd_WKs, mean_TD_ind,
  sd_TD_ind, mean_REVtd, sd_REVtd) %>%
  summarise(
    DAYs_Index_mean     = mean(DAYs_index, na.rm = TRUE),
    DAYs_Index_sd     = sd(DAYs_index, na.rm = TRUE),
    WKs_Index_mean     = mean(WKs_index, na.rm = TRUE),
    WKs_Index_sd     = sd(WKs_index, na.rm = TRUE),
    REVtd_Index_mean = mean(REVtd_index, na.rm = TRUE),
    REVtd_Index_sd = sd(REVtd_index, na.rm = TRUE),
    .groups = "drop"
  )

# Left-join group averages into sig_pairs_mean_y_neat
sig_pairs_with_group_means <- halves_neat %>%
  left_join(DatMean,
            by = c("SPECIES", "PROJECT"))

## Estimate Coefficient of Variation–like Metric
#Normalizes by standard deviation of all values.
#Interpreted as “how many SDs apart are these two values.”
unique(sig_pairs_with_group_means$y_var)

# Define a lookup between y_var values and their corresponding sd columns
sd_lookup <- c(
  DAYs_index  = "DAYs_Index_sd",
  WKs_index    = "WKs_Index_sd",
  REVtd_index    = "REVtd_Index_sd"
)

sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  rowwise() %>%
  mutate(
    CVdiff = {
      sd_col <- sd_lookup[[y_var]]
      if (is.null(sd_col)) NA_real_
      else abs(Diff_range) / get(sd_col)
    }
  ) %>%
  ungroup()


sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  mutate(CVdiff = round(CVdiff,3))

colnames(sig_pairs_with_group_means)
sig_pairs_with_group_means <- sig_pairs_with_group_means %>%
  select(-c(17:30))

# Save and pretty-print
write_csv(sig_pairs_with_group_means, file.path(out_dir, "TD_significant_pairs_mean_y_by_halves.csv"))

# ==================== WORD REPORT (OFFICER) =======================================
docx_path <- file.path(out_dir, "TD_Correlations_Report.docx")

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

doc <- body_add_par(doc, "TD — Correlations & Half-Period Comparison", style = "heading 1")
doc <- body_add_par(doc, "Significant pairs & first vs second period comparison", style = "heading 2")
doc <- body_add_table(doc, value = Comparisons, style = "Normal Table")

# Add the “all groups” pair-plot files if you generated them
doc <- body_add_par(doc, "Significant pairs — facet plots (★ marks significant panels)", style = "heading 2")
for (img in list.files(out_dir, pattern = "^TD_.*\\.png$", full.names = TRUE)) {
  doc <- body_add_img(doc, src = img, width = 7, height = 4.7)
}

print(doc, target = docx_path)
message("Word report saved at: ", normalizePath(docx_path, winslash = "/", mustWork = FALSE))

