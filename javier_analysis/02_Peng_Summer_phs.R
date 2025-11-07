##### Review Penguin indices
##### generate SUMMER INDICES
##### post-hatch success (phs)

## Watters Code:
# post-hatch success (phs) (numbers of chicks creched/numbers of chicks hatched)
# bigger indicates better summer
phs <- read.csv("./Supplementary Files/success.csv") # post-hatch success
phs$phs<-phs$N_CRECHE/phs$N_CHICKS  # No. creche chicks per chick hatched
phs$phs<-log(phs$phs/(1-phs$phs))  # standardized to 0
phs$matchme<-paste(phs$PROJECT,phs$SPECIES,sep="|")
tt<-tapply(phs$phs,list(phs$matchme),mean,na.rm=TRUE) # mean per site and species
ttt<-tapply(phs$phs,list(phs$matchme),sd,na.rm=TRUE)
mean.phs<-tt[match(phs$matchme,names(tt))]
sd.phs<-ttt[match(phs$matchme,names(ttt))]
phs$std.logit.phs<-(phs$phs-mean.phs)/sd.phs  # INDEX
phs<-phs[,-c(4:9)]
names(phs)[4]<-"index"
phs$param=rep("PHS",dim(phs)[1])
phs$season=rep("S",dim(phs)[1])
# summer indices are relevant to the second year in the split-season designation
phs$cal.yr<-as.numeric(substr(phs$YEAR,1,4))+1

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
stopifnot(exists("phs"), is.data.frame(phs))

# --- Output locations ---------------------------
out_dir <- "./javier_analysis/review_PHS"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Data Prep ----------------------------------------
# Read the same source path used by the Watters code
ade <- read.csv("./Supplementary Files/success.csv", header=TRUE, stringsAsFactors = FALSE) # post-hatch success

ade <- ade |>
  mutate(
    # DATE        = suppressWarnings(lubridate::mdy(DATE)),
    cal.yr      = as.integer(substr(YEAR, 1, 4))+1,
    # julian_days = lubridate::yday(DATE)
    PHS   = (N_CRECHE/N_CHICKS),  # post-hatching success
    PHS_ind = log(PHS / (1 - PHS)),  # standardized to 0
    BS    = (N_CRECHE/N_EGGS),   # breeding success
    BS_ind = log(BS / (1 - BS))  # standardized to 0
  )

# Step 1: YEAR x PROJECT x SPECIES means
step1 <- ade %>%
  filter(!is.na(N_CRECHE)) %>%
  group_by(cal.yr, PROJECT, SPECIES) %>%
  summarise(
    PHSmean   = mean(PHS, na.rm = TRUE),
    BSmean   = mean(BS, na.rm = TRUE),
    PHS_ind_mean  = mean(PHS_ind, na.rm = TRUE),
    BS_ind_mean  = mean(BS_ind, na.rm = TRUE),
    .groups  = "drop"
  )

# Step 2: Within PROJECT x SPECIES: global mean/sd and day mean
step2 <- step1 %>%
  group_by(PROJECT, SPECIES) %>%
  summarise(
    mean_PHS = mean(PHSmean, na.rm = TRUE),
    sd_PHS   = sd(PHSmean, na.rm = TRUE),
    mean_BS = mean(BSmean, na.rm = TRUE),
    sd_BS   = sd(BSmean, na.rm = TRUE),
    mean_PHS_ind = mean(PHS_ind_mean, na.rm = TRUE),
    sd_PHS_ind   = sd(PHS_ind_mean, na.rm = TRUE),
    mean_BS_ind = mean(BS_ind_mean, na.rm = TRUE),
    sd_BS_ind   = sd(BS_ind_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Join & indices (guard against sd_WT == 0)
ADE_sum <- step1 %>%
  left_join(step2, by = c("PROJECT","SPECIES")) %>%
  mutate(
    PH_index     = ifelse(sd_PHS_ind > 0, (PHS_ind_mean - mean_PHS_ind)/sd_PHS_ind, NA_real_),
    BS_index = ifelse(sd_BS_ind > 0, (BS_ind_mean - mean_BS_ind)/sd_BS_ind, NA_real_),
  ) %>%
  arrange(SPECIES, PROJECT, cal.yr)

# -------------------- CORRELATIONS ------------------------------
# Robust correlation-by-group helper with status flag
source(file.path("./javier_analysis/corr_tbl_ts.R"))

# choose candidate variables (means)
colnames(ADE_sum)
cand <- c("cal.yr", "PHSmean",  "BSmean", "PHS_ind_mean", "BS_ind_mean",
          "PH_index", "BS_index")

cand <- intersect(cand, names(ADE_sum))
pairs <- combn(cand, 2, simplify = FALSE)

run_corr_pair <- function(df, nm_x, nm_y) {
  x_sym <- sym(nm_x); y_sym <- sym(nm_y)
  corr_tbl_ts(df, !!x_sym, !!y_sym) %>% mutate(x_var=nm_x, y_var=nm_y, .before=1)
}

all_corrs   <- map_dfr(pairs, ~ run_corr_pair(ADE_sum,   .x[1], .x[2]))

# Save
write_csv(all_corrs, file.path(out_dir, "PHS_all_combinations_corr.csv"))

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
    message("No significant pairs found in PHS")
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
    
    fpath <- file.path(out_dir, paste0("PHS_", sanitize(y_var), "_", sanitize(x_var), ".png"))
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
    message("No pairs found in PHS")
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
    
    fpath <- file.path(out_dir, paste0("PHS_", sanitize(y_var), "_", sanitize(x_var), ".png"))
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
DatMean <- ADE_sum %>%
  group_by(SPECIES, PROJECT, mean_PHS, sd_PHS, mean_BS, sd_BS, mean_BS_ind, sd_BS_ind, mean_PHS_ind, sd_PHS_ind) %>%
  summarise(
    PH_Index_mean     = mean(PH_index, na.rm = TRUE),
    PH_Index_sd     = sd(PH_index, na.rm = TRUE),
    BS_Index_mean = mean(BS_index, na.rm = TRUE),
    BS_Index_sd = sd(BS_index, na.rm = TRUE),
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
  PH_index  = "PH_Index_sd",
  PHS_ind_mean    = "sd_PHS_ind",
  BS_ind_mean    = "sd_BS_ind",
  BS_index     = "BS_Index_sd",
  BSmean  = "sd_BS"
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
  select(-mean_PHS, -sd_PHS, -mean_BS, -sd_BS, -mean_BS_ind, -sd_BS_ind, -mean_PHS_ind, -sd_PHS_ind,   
         -PH_Index_mean, -PH_Index_sd, -BS_Index_mean, -BS_Index_sd)

# Save and pretty-print
write_csv(sig_pairs_with_group_means, file.path(out_dir, "PHS_significant_pairs_mean_y_by_halves.csv"))

# ==================== WORD REPORT (OFFICER) =======================================
docx_path <- file.path(out_dir, "PHS_Correlations_Report.docx")

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

doc <- body_add_par(doc, "PHS — Correlations & Half-Period Comparison", style = "heading 1")
doc <- body_add_par(doc, "Significant pairs & first vs second period comparison", style = "heading 2")
doc <- body_add_table(doc, value = Comparisons, style = "Normal Table")

# Add the “all groups” pair-plot files if you generated them
doc <- body_add_par(doc, "Significant pairs — facet plots (★ marks significant panels)", style = "heading 2")
for (img in list.files(out_dir, pattern = "^PHS_.*\\.png$", full.names = TRUE)) {
  doc <- body_add_img(doc, src = img, width = 7, height = 4.7)
}

print(doc, target = docx_path)
message("Word report saved at: ", normalizePath(docx_path, winslash = "/", mustWork = FALSE))

