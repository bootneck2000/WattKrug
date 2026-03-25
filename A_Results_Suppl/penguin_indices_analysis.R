# =============================================================================
# Penguin Indices Analysis — Model01_filtered
# Data: A_Results_Suppl/Model01_filtered/junk.rds
# Analyses:
#   1. Boxplots: index ~ sam.sign | oni.class | bclass | hrclass (by season)
#   2. Statistical tests: Kruskal-Wallis / Mann-Whitney U per category
#   3. Scatter plots: index ~ survey (LKB) and LHR (by season)
#   4. Spearman correlations: index vs LKB and LHR
#   5. Summary results table (N, p-values)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)

# ------------------------------------------------------------------------------
# 0. Load & prepare data
# ------------------------------------------------------------------------------
dat <- readRDS("./A_Results_Suppl/Model01_filtered/junk.rds")

# Compute LHR (harvest rate = catch / LKB survey biomass)
dat <- dat %>%
  mutate(
    LHR     = catch / survey,
    # Convert categorical variables to ordered factors
    sam.sign  = factor(sam.sign,  levels = c("Neg", "Pos")),
    oni.class = factor(oni.class, levels = c("Cool", "Neutral", "Warm")),
    bclass    = factor(bclass,    levels = c(1, 2),
                       labels = c("Low", "High")),
    hrclass   = factor(hrclass,   levels = c(1, 2, 3),
                       labels = c("Low", "Mid", "High")),
    season    = factor(season,    levels = c("S", "W"),
                       labels = c("Summer", "Winter"))
  )

# Season-specific param sets
summer_params <- c("FWT", "PHS", "REVTD")
winter_params <- c("MML", "FML", "EGG", "REVCID", "REC")

dat_summer <- dat %>% filter(season == "Summer", param %in% summer_params)
dat_winter <- dat %>% filter(season == "Winter", param %in% winter_params)

cat("Summer rows:", nrow(dat_summer), "| Winter rows:", nrow(dat_winter), "\n")
cat("Summer params present:", paste(unique(dat_summer$param), collapse = ", "), "\n")
cat("Winter params present:", paste(unique(dat_winter$param), collapse = ", "), "\n")

# Colour palette per param
param_colors <- c(
  FWT    = "#2196F3", PHS    = "#4CAF50", REVTD  = "#FF9800",
  MML    = "#9C27B0", FML    = "#F44336", EGG    = "#00BCD4",
  REVCID = "#795548", REC    = "#607D8B"
)

# Helper: strip + theme consistent across plots
theme_ark <- function() {
  theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#E3F2FD"),
      strip.text       = element_text(face = "bold", size = 10),
      axis.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# 1. BOXPLOTS — index vs categorical variables
# =============================================================================

cat_vars    <- c("sam.sign", "oni.class", "bclass", "hrclass")
cat_labels  <- c("SAM Sign", "ONI Class", "Biomass Class", "Harvest Rate Class")

make_boxplot <- function(df, cat_var, cat_label, season_label) {
  ggplot(df, aes(x = .data[[cat_var]], y = index,
                 fill = param, colour = param)) +
    geom_boxplot(alpha = 0.55, outlier.shape = 21, outlier.size = 1.8,
                 position = position_dodge(width = 0.8), width = 0.6) +
    scale_fill_manual(values   = param_colors, name = "Index") +
    scale_colour_manual(values = param_colors, name = "Index") +
    labs(
      title = paste0(season_label, " — Index vs ", cat_label),
      x     = cat_label,
      y     = "Penguin Index"
    ) +
    theme_ark()
}

# --- Summer boxplots ---
bp_summer <- lapply(seq_along(cat_vars), function(i) {
  make_boxplot(dat_summer, cat_vars[i], cat_labels[i], "Summer")
})
names(bp_summer) <- cat_vars

fig_bp_summer <- wrap_plots(bp_summer, ncol = 2) +
  plot_annotation(
    title    = "Summer Penguin Indices vs Categorical Drivers",
    subtitle = "Params: FWT, PHS, REVTD",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 11))
  ); fig_bp_summer

ggsave("./A_Results_Suppl/Model01_filtered/Fig_Boxplot_Summer.png",
       fig_bp_summer, width = 14, height = 10, dpi = 180)
cat("Saved: Fig_Boxplot_Summer.png\n")

# --- Winter boxplots ---
if (nrow(dat_winter) > 0) {
  bp_winter <- lapply(seq_along(cat_vars), function(i) {
    make_boxplot(dat_winter, cat_vars[i], cat_labels[i], "Winter")
  })
  names(bp_winter) <- cat_vars

  fig_bp_winter <- wrap_plots(bp_winter, ncol = 2) +
    plot_annotation(
      title    = "Winter Penguin Indices vs Categorical Drivers",
      subtitle = "Params: MML, FML, EGG, REVCID (REC if available)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                       plot.subtitle = element_text(size = 11))
    )

  ggsave("A_Results_Suppl/Model01_filtered/Fig_Boxplot_Winter.png",
         fig_bp_winter, width = 14, height = 10, dpi = 180)
  cat("Saved: Fig_Boxplot_Winter.png\n")
}
fig_bp_winter
# =============================================================================
# 2. STATISTICAL TESTS — index vs categorical variables
# =============================================================================
# For 2-level vars (sam.sign, bclass): Wilcoxon rank-sum test
# For 3-level vars (oni.class, hrclass): Kruskal-Wallis + Dunn post-hoc

run_cat_tests <- function(df, season_label) {

  results <- list()

  for (cat_var in cat_vars) {
    n_levels <- nlevels(df[[cat_var]])

    for (p in unique(df$param)) {
      sub <- df %>% filter(param == p) %>%
        select(index, !!sym(cat_var)) %>%
        drop_na()

      if (nrow(sub) < 4) next

      grp_n <- sub %>% group_by(.data[[cat_var]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(label = paste0(.data[[cat_var]], " (n=", n, ")")) %>%
        pull(label) %>% paste(collapse = "; ")

      if (n_levels == 2) {
        # Wilcoxon rank-sum (Mann-Whitney U)
        grps  <- levels(df[[cat_var]])
        x1    <- sub %>% filter(.data[[cat_var]] == grps[1]) %>% pull(index)
        x2    <- sub %>% filter(.data[[cat_var]] == grps[2]) %>% pull(index)
        if (length(x1) < 2 || length(x2) < 2) next
        wt    <- wilcox.test(x1, x2, exact = FALSE)
        results[[length(results) + 1]] <- tibble(
          season    = season_label,
          param     = p,
          predictor = cat_var,
          test      = "Wilcoxon",
          comparison = paste(grps, collapse = " vs "),
          N          = nrow(sub),
          groups_N   = grp_n,
          statistic  = round(wt$statistic, 3),
          p_value    = round(wt$p.value, 4),
          sig        = ifelse(wt$p.value < 0.001, "***",
                       ifelse(wt$p.value < 0.01,  "**",
                       ifelse(wt$p.value < 0.05,  "*",
                       ifelse(wt$p.value < 0.10,  ".", "ns"))))
        )

      } else {
        # Kruskal-Wallis overall
        kw <- kruskal.test(reformulate(cat_var, response = "index"), data = sub)
        results[[length(results) + 1]] <- tibble(
          season     = season_label,
          param      = p,
          predictor  = cat_var,
          test       = "Kruskal-Wallis",
          comparison = "Overall",
          N          = nrow(sub),
          groups_N   = grp_n,
          statistic  = round(kw$statistic, 3),
          p_value    = round(kw$p.value, 4),
          sig        = ifelse(kw$p.value < 0.001, "***",
                       ifelse(kw$p.value < 0.01,  "**",
                       ifelse(kw$p.value < 0.05,  "*",
                       ifelse(kw$p.value < 0.10,  ".", "ns"))))
        )

        # Dunn post-hoc pairwise (if K-W significant or always for completeness)
        if (kw$p.value < 0.10 && length(unique(sub[[cat_var]])) > 1) {
          dunn <- dunn_test(sub, as.formula(paste("index ~", cat_var)),
                            p.adjust.method = "BH")
          for (k in seq_len(nrow(dunn))) {
            results[[length(results) + 1]] <- tibble(
              season     = season_label,
              param      = p,
              predictor  = cat_var,
              test       = "Dunn (BH)",
              comparison = paste(dunn$group1[k], "vs", dunn$group2[k]),
              N          = nrow(sub),
              groups_N   = grp_n,
              statistic  = round(dunn$statistic[k], 3),
              p_value    = round(dunn$p.adj[k], 4),
              sig        = ifelse(dunn$p.adj[k] < 0.001, "***",
                           ifelse(dunn$p.adj[k] < 0.01,  "**",
                           ifelse(dunn$p.adj[k] < 0.05,  "*",
                           ifelse(dunn$p.adj[k] < 0.10,  ".", "ns"))))
            )
          }
        }
      }
    }
  }

  bind_rows(results)
}

stats_summer <- run_cat_tests(dat_summer, "Summer")
stats_winter <- if (nrow(dat_winter) > 0) run_cat_tests(dat_winter, "Winter") else tibble()
stats_cat    <- bind_rows(stats_summer, stats_winter)

cat("\n--- Categorical test results (significant only, p < 0.10) ---\n")
print(stats_cat %>% filter(sig != "ns") %>%
        select(season, param, predictor, test, comparison, N, p_value, sig), n = 40)

# =============================================================================
# 3. SCATTER PLOTS — index vs LKB (survey) and LHR
# =============================================================================

make_scatter <- function(df, x_var, x_label, season_label) {
  ggplot(df, aes(x = .data[[x_var]], y = index,
                 colour = param, fill = param)) +
    geom_point(alpha = 0.65, size = 2.2, shape = 21, stroke = 0.4) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
    facet_wrap(~param, scales = "free_x", ncol = 3) +
    scale_colour_manual(values = param_colors, guide = "none") +
    scale_fill_manual(values   = param_colors, guide = "none") +
    labs(
      title    = paste0(season_label, " — Index vs ", x_label),
      x        = x_label,
      y        = "Penguin Index"
    ) +
    theme_ark() +
    theme(strip.text = element_text(face = "bold"))
}

# Summer LKB
fig_sc_summer_lkb <- make_scatter(dat_summer, "survey", "Krill Biomass (LKB, tonnes)", "Summer")
ggsave("A_Results_Suppl/Model01_filtered/Fig_Scatter_Summer_LKB.png",
       fig_sc_summer_lkb, width = 12, height = 5, dpi = 180)
cat("Saved: Fig_Scatter_Summer_LKB.png\n")

# Summer LHR
fig_sc_summer_lhr <- make_scatter(dat_summer, "LHR", "Harvest Rate (LHR = catch/LKB)", "Summer")
ggsave("A_Results_Suppl/Model01_filtered/Fig_Scatter_Summer_LHR.png",
       fig_sc_summer_lhr, width = 12, height = 5, dpi = 180)
cat("Saved: Fig_Scatter_Summer_LHR.png\n")

# Winter LKB & LHR (if data available)
if (nrow(dat_winter) > 0) {
  fig_sc_winter_lkb <- make_scatter(dat_winter, "survey", "Krill Biomass (LKB, tonnes)", "Winter")
  ggsave("A_Results_Suppl/Model01_filtered/Fig_Scatter_Winter_LKB.png",
         fig_sc_winter_lkb, width = 10, height = 5, dpi = 180)
  cat("Saved: Fig_Scatter_Winter_LKB.png\n")

  fig_sc_winter_lhr <- make_scatter(dat_winter, "LHR", "Harvest Rate (LHR = catch/LKB)", "Winter")
  ggsave("./A_Results_Suppl/Model01_filtered/Fig_Scatter_Winter_LHR.png",
         fig_sc_winter_lhr, width = 10, height = 5, dpi = 180)
  cat("Saved: Fig_Scatter_Winter_LHR.png\n")
}

# =============================================================================
# 4. SPEARMAN CORRELATIONS — index vs LKB and LHR
# =============================================================================

run_corr_tests <- function(df, season_label) {
  cont_vars  <- c("survey", "LHR")
  cont_labels <- c("LKB (survey biomass)", "LHR (catch/LKB)")
  results <- list()

  for (i in seq_along(cont_vars)) {
    for (p in unique(df$param)) {
      sub <- df %>% filter(param == p) %>%
        select(index, val = !!sym(cont_vars[i])) %>%
        drop_na() %>%
        filter(is.finite(val), is.finite(index))

      if (nrow(sub) < 5) next

      ct <- cor.test(sub$index, sub$val,
                     method = "spearman", exact = FALSE)
      results[[length(results) + 1]] <- tibble(
        season     = season_label,
        param      = p,
        predictor  = cont_vars[i],
        predictor_label = cont_labels[i],
        N          = nrow(sub),
        rho        = round(ct$estimate, 3),
        p_value    = round(ct$p.value, 4),
        sig        = ifelse(ct$p.value < 0.001, "***",
                     ifelse(ct$p.value < 0.01,  "**",
                     ifelse(ct$p.value < 0.05,  "*",
                     ifelse(ct$p.value < 0.10,  ".", "ns"))))
      )
    }
  }
  bind_rows(results)
}

corr_summer <- run_corr_tests(dat_summer, "Summer")
corr_winter <- if (nrow(dat_winter) > 0) run_corr_tests(dat_winter, "Winter") else tibble()
corr_all    <- bind_rows(corr_summer, corr_winter)

cat("\n--- Spearman correlations ---\n")
print(corr_all %>% select(season, param, predictor, N, rho, p_value, sig), n = 40)

# =============================================================================
# 5. SUMMARY RESULTS TABLE (combined)
# =============================================================================

# — Categorical tests table (clean format)
table_cat <- stats_cat %>%
  select(Season = season,
         Param  = param,
         Predictor = predictor,
         Test   = test,
         Comparison = comparison,
         N,
         `Groups (N)` = groups_N,
         Statistic = statistic,
         `p-value`  = p_value,
         Sig = sig) %>%
  arrange(Season, Predictor, Param, Test)

# — Continuous correlations table
table_corr <- corr_all %>%
  mutate(Test       = "Spearman rho",
         Comparison = NA_character_,
         `Groups (N)` = as.character(N)) %>%
  select(Season      = season,
         Param       = param,
         Predictor   = predictor_label,
         Test,
         Comparison,
         N,
         `Groups (N)`,
         Statistic   = rho,
         `p-value`   = p_value,
         Sig         = sig) %>%
  arrange(Season, Predictor, Param)

table_full <- bind_rows(table_cat, table_corr)

# Save as CSV
write_csv(table_cat,  "./A_Results_Suppl/Model01_filtered/Table_CatTests.csv")
write_csv(table_corr, "./A_Results_Suppl/Model01_filtered/Table_CorrTests.csv")
write_csv(table_full, "./A_Results_Suppl/Model01_filtered/Table_AllResults.csv")

cat("\n--- Full results table saved ---\n")
cat("  Table_CatTests.csv  — categorical test results\n")
cat("  Table_CorrTests.csv — Spearman correlation results\n")
cat("  Table_AllResults.csv — combined\n")

# Print summaries to console
cat("\n=== CATEGORICAL TESTS SUMMARY ===\n")
print(table_cat, n = 80)
cat("\n=== CORRELATION TESTS SUMMARY ===\n")
print(table_corr, n = 40)

cat("\n\nAll done. Output files in: A_Results_Suppl/Model01_filtered/\n")
