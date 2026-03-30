#### RELATIVE IMPORTANCE OF ONI AND LHR ON PENGUIN PERFORMANCE
## Uses saved MCMC posteriors from Model 02 (SSMUs) — no model re-run needed.
##
## Model recap (sum-to-zero contrasts):
##   beta[1] = intercept
##   beta[2] = LKB class effect  (low vs high, >1 Mt)
##   beta[3] = LHR medium effect (0.01 < LHR < 0.1)
##   beta[4] = LHR high effect   (LHR >= 0.1)
##   beta[5] = ONI neutral       (-0.5 <= ONI < 0.5)
##   beta[6] = ONI warm          (ONI >= 0.5)
##
## Sum-to-zero level effects:
##   ONI: cool = -beta5 - beta6;  neutral = +beta5;  warm = +beta6
##   LHR: low  = -beta3 - beta4;  med     = +beta3;  high = +beta4
##   LKB: low  = -beta2;          high    = +beta2

library(tidyverse)
library(coda)
library(ggmcmc)

out.dir <- "./A_Results/Model02_SSMU/"
res.dir <- "./A_Results/Relative_Importance/"
dir.create(res.dir, showWarnings = FALSE, recursive = TRUE)

## ── Load saved MCMC results ────────────────────────────────────────────────
hr.params.post  <- readRDS(file.path(out.dir, "hr.params.post.rds"))
hr.derived.post <- readRDS(file.path(out.dir, "hr.derived.post.rds"))

hr.params.s  <- ggs(hr.params.post)
hr.derived.s <- ggs(hr.derived.post)

## ── 1.  Beta posteriors → level effects ───────────────────────────────────
beta_wide <- hr.params.s %>%
  filter(grepl("^beta\\[", as.character(Parameter))) %>%
  mutate(beta_idx = as.integer(gsub("beta\\[|\\]", "", as.character(Parameter)))) %>%
  select(Iteration, Chain, beta_idx, value) %>%
  pivot_wider(names_from = beta_idx, values_from = value, names_prefix = "beta") %>%
  mutate(
    # ONI level deviations (sum-to-zero)
    oni_cool    = -beta5 - beta6,
    oni_neutral =  beta5,
    oni_warm    =  beta6,
    # LHR level deviations (sum-to-zero)
    lhr_low     = -beta3 - beta4,
    lhr_med     =  beta3,
    lhr_high    =  beta4,
    # LKB level deviations (sum-to-zero)
    lkb_low     = -beta2,
    lkb_high    =  beta2,
    # Effect range per MCMC sample (max - min across levels)
    oni_range = pmax(oni_cool, oni_neutral, oni_warm) - pmin(oni_cool, oni_neutral, oni_warm),
    lhr_range = pmax(lhr_low,  lhr_med,    lhr_high)  - pmin(lhr_low,  lhr_med,    lhr_high),
    lkb_range = abs(2 * beta2)   # two-level factor: range = 2|beta2|
  )

## ── 2.  Summary statistics ─────────────────────────────────────────────────
range_summary <- beta_wide %>%
  select(oni_range, lhr_range, lkb_range) %>%
  pivot_longer(everything(), names_to = "Predictor", values_to = "Range") %>%
  mutate(Predictor = recode(Predictor,
    oni_range = "ONI",
    lhr_range = "LHR",
    lkb_range = "LKB (krill biomass)"
  )) %>%
  group_by(Predictor) %>%
  summarise(
    Median = round(median(Range), 3),
    Mean   = round(mean(Range),   3),
    SD     = round(sd(Range),     3),
    Q2.5   = round(quantile(Range, 0.025), 3),
    Q97.5  = round(quantile(Range, 0.975), 3),
    .groups = "drop"
  )

cat("\n=== Effect-range summary (standardised performance units) ===\n")
print(range_summary)

p_lhr_gt_oni <- mean(beta_wide$lhr_range > beta_wide$oni_range)
p_oni_gt_lkb <- mean(beta_wide$oni_range > beta_wide$lkb_range)
p_lhr_gt_lkb <- mean(beta_wide$lhr_range > beta_wide$lkb_range)

cat(sprintf("\nP(LHR range > ONI range) = %.3f\n", p_lhr_gt_oni))
cat(sprintf("P(ONI range > LKB range) = %.3f\n", p_oni_gt_lkb))
cat(sprintf("P(LHR range > LKB range) = %.3f\n", p_lhr_gt_lkb))

## ── 3.  Model-computed comparative probabilities (prob[7:12]) ──────────────
# prob[7]:  P(mu | med LHR  < mu | neutral ONI)  — med LHR worse than neutral ONI
# prob[8]:  P(mu | high LHR < mu | neutral ONI)  — high LHR worse than neutral ONI
# prob[9]:  P(mu | high LKB < mu | neutral ONI)
# prob[10]: P(mu | med LHR  < mu | warm ONI)     — med LHR worse than warm ONI
# prob[11]: P(mu | high LHR < mu | warm ONI)     — high LHR worse than warm ONI
# prob[12]: P(mu | high LKB < mu | warm ONI)

prob_labels <- c(
  "prob[1]"  = "P(high LKB reduces perf. vs reference)",
  "prob[2]"  = "P(med LHR reduces perf. vs reference)",
  "prob[3]"  = "P(high LHR reduces perf. vs reference)",
  "prob[4]"  = "P(neutral ONI reduces perf. vs reference)",
  "prob[5]"  = "P(warm ONI reduces perf. vs reference)",
  "prob[6]"  = "P(worst case reduces perf. vs reference)",
  "prob[7]"  = "P(med LHR effect worse than neutral ONI)",
  "prob[8]"  = "P(high LHR effect worse than neutral ONI)",
  "prob[9]"  = "P(high LKB effect worse than neutral ONI)",
  "prob[10]" = "P(med LHR effect worse than warm ONI)",
  "prob[11]" = "P(high LHR effect worse than warm ONI)",
  "prob[12]" = "P(high LKB effect worse than warm ONI)"
)

prob_summ <- hr.derived.s %>%
  filter(as.character(Parameter) %in% names(prob_labels)) %>%
  group_by(Parameter) %>%
  summarise(Probability = round(mean(value), 3), .groups = "drop") %>%
  mutate(Label = prob_labels[as.character(Parameter)]) %>%
  arrange(Parameter)

cat("\n=== Comparative probabilities (from JAGS model) ===\n")
print(prob_summ %>% select(Parameter, Label, Probability), n = Inf)

write.csv(range_summary, file.path(res.dir, "effect_range_summary.csv"),    row.names = FALSE)
write.csv(prob_summ,     file.path(res.dir, "comparative_probabilities.csv"), row.names = FALSE)

## ── 4.  Figure 1: Posterior densities of effect ranges ────────────────────
range_long <- beta_wide %>%
  select(oni_range, lhr_range, lkb_range) %>%
  pivot_longer(everything(), names_to = "Predictor", values_to = "Range") %>%
  mutate(Predictor = recode(Predictor,
    oni_range = "ONI",
    lhr_range = "LHR",
    lkb_range = "LKB (krill biomass)"
  ))

fig1 <- ggplot(range_long, aes(x = Range, fill = Predictor, colour = Predictor)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = c("ONI" = "#2171b5", "LHR" = "#cb181d", "LKB (krill biomass)" = "#238b45")) +
  scale_colour_manual(values = c("ONI" = "#2171b5", "LHR" = "#cb181d", "LKB (krill biomass)" = "#238b45")) +
  labs(
    title    = "Relative importance: effect ranges on penguin performance",
    subtitle = paste0("P(LHR range > ONI range) = ", round(p_lhr_gt_oni, 3),
                      "   |   P(ONI range > LKB range) = ", round(p_oni_gt_lkb, 3)),
    x      = "Effect range (standardised performance units)",
    y      = "Posterior density",
    fill   = "Predictor",
    colour = "Predictor"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(res.dir, "Fig1_effect_ranges.png"), fig1, width = 8, height = 5, dpi = 300)

## ── 5.  Figure 2: Level effects for ONI and LHR ───────────────────────────
level_long <- beta_wide %>%
  select(oni_cool, oni_neutral, oni_warm, lhr_low, lhr_med, lhr_high) %>%
  pivot_longer(everything(), names_to = "Level", values_to = "Effect") %>%
  mutate(
    Predictor = if_else(str_starts(Level, "oni"), "ONI", "LHR"),
    Level = recode(Level,
      oni_cool    = "Cool (ONI <= -0.5)",
      oni_neutral = "Neutral (-0.5 < ONI < 0.5)",
      oni_warm    = "Warm (ONI >= 0.5)",
      lhr_low     = "Low  (LHR <= 0.01)",
      lhr_med     = "Medium (0.01 < LHR < 0.1)",
      lhr_high    = "High (LHR >= 0.1)"
    ),
    Level = factor(Level, levels = c(
      "Cool (ONI <= -0.5)", "Neutral (-0.5 < ONI < 0.5)", "Warm (ONI >= 0.5)",
      "Low  (LHR <= 0.01)", "Medium (0.01 < LHR < 0.1)", "High (LHR >= 0.1)"
    ))
  )

fig2 <- ggplot(level_long, aes(x = Effect, fill = Level, colour = Level)) +
  geom_density(alpha = 0.40, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ Predictor, ncol = 1, scales = "free_y") +
  labs(
    title    = "Posterior distributions of predictor-level effects",
    subtitle = "Level effects derived from sum-to-zero beta posteriors",
    x      = "Effect on performance (standardised units)",
    y      = "Posterior density",
    fill   = "Level",
    colour = "Level"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(res.dir, "Fig2_level_effects.png"), fig2, width = 9, height = 7, dpi = 300)

## ── 6.  Figure 3: Marginal effects via mu.index.new ───────────────────────
# pred.matrix row mapping (from Mod_02_SSMUs.R):
#  1: cool,    low LKB, low LHR  ← ONI reference for LHR comparison
#  3: cool,    low LKB, med LHR
#  5: cool,    low LKB, high LHR
#  7: neutral, low LKB, low LHR  ← LHR reference for ONI comparison
# 13: warm,    low LKB, low LHR

mu_s <- hr.derived.s %>%
  filter(grepl("^mu\\.index\\.new\\[", as.character(Parameter))) %>%
  mutate(row_idx = as.integer(gsub("mu\\.index\\.new\\[|\\]", "", as.character(Parameter))))

oni_marg <- mu_s %>%
  filter(row_idx %in% c(1, 7, 13)) %>%
  mutate(ONI = recode(as.character(row_idx),
    "1"  = "Cool (ONI <= -0.5)",
    "7"  = "Neutral",
    "13" = "Warm (ONI >= 0.5)"
  ))

lhr_marg <- mu_s %>%
  filter(row_idx %in% c(1, 3, 5)) %>%
  mutate(LHR = recode(as.character(row_idx),
    "1" = "Low (LHR <= 0.01)",
    "3" = "Medium (0.01 < LHR < 0.1)",
    "5" = "High (LHR >= 0.1)"
  ))

fig3a <- ggplot(oni_marg, aes(x = value, fill = ONI, colour = ONI)) +
  geom_density(alpha = 0.40, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = c("Cool (ONI <= -0.5)" = "#2171b5",
                                  "Neutral"            = "#6baed6",
                                  "Warm (ONI >= 0.5)"  = "#fd8d3c")) +
  scale_colour_manual(values = c("Cool (ONI <= -0.5)" = "#2171b5",
                                  "Neutral"            = "#6baed6",
                                  "Warm (ONI >= 0.5)"  = "#fd8d3c")) +
  labs(title    = "ONI marginal effect  [LKB = low, LHR = low]",
       x = "Expected performance", y = "Posterior density",
       fill = "ONI class", colour = "ONI class") +
  theme_bw(base_size = 12)

fig3b <- ggplot(lhr_marg, aes(x = value, fill = LHR, colour = LHR)) +
  geom_density(alpha = 0.40, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = c("Low (LHR <= 0.01)"          = "#74c476",
                                  "Medium (0.01 < LHR < 0.1)" = "#fd8d3c",
                                  "High (LHR >= 0.1)"          = "#cb181d")) +
  scale_colour_manual(values = c("Low (LHR <= 0.01)"          = "#74c476",
                                  "Medium (0.01 < LHR < 0.1)" = "#fd8d3c",
                                  "High (LHR >= 0.1)"          = "#cb181d")) +
  labs(title    = "LHR marginal effect  [ONI = cool, LKB = low]",
       x = "Expected performance", y = "Posterior density",
       fill = "LHR class", colour = "LHR class") +
  theme_bw(base_size = 12)

ggsave(file.path(res.dir, "Fig3a_ONI_marginal.png"), fig3a, width = 8, height = 4.5, dpi = 300)
ggsave(file.path(res.dir, "Fig3b_LHR_marginal.png"), fig3b, width = 8, height = 4.5, dpi = 300)

## ── 7.  Figure 4: Comparative probability bar chart ───────────────────────
# Highlight ONI vs LHR comparisons (prob[7:12])
prob_plot <- prob_summ %>%
  mutate(Group = case_when(
    as.character(Parameter) %in% c("prob[7]","prob[8]","prob[10]","prob[11]") ~ "LHR vs ONI",
    as.character(Parameter) %in% c("prob[9]","prob[12]")                      ~ "LKB vs ONI",
    TRUE                                                                       ~ "vs reference"
  ))

fig4 <- ggplot(prob_plot, aes(x = Probability, y = reorder(Label, Probability), fill = Group)) +
  geom_col(alpha = 0.75) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = c("LHR vs ONI"   = "#cb181d",
                                "LKB vs ONI"   = "#238b45",
                                "vs reference" = "#969696")) +
  labs(
    title    = "Comparative probabilities: LHR and ONI effects",
    subtitle = "P(effect of row < effect of column) from posterior expectations",
    x    = "Posterior probability",
    y    = NULL,
    fill = "Comparison"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(res.dir, "Fig4_comparative_probs.png"), fig4, width = 10, height = 5.5, dpi = 300)

## ── 8.  Figure 5: LHR range vs ONI range scatter (posterior samples) ──────
fig5 <- ggplot(beta_wide %>% sample_n(min(3000, n())),
               aes(x = oni_range, y = lhr_range)) +
  geom_point(alpha = 0.08, size = 0.5, colour = "#2171b5") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = 0.8) +
  labs(
    title    = "LHR range vs ONI range  (posterior samples)",
    subtitle = paste0("Points above dashed line: LHR has larger effect range than ONI\n",
                      "P(LHR range > ONI range) = ", round(p_lhr_gt_oni, 3)),
    x = "ONI effect range",
    y = "LHR effect range"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(res.dir, "Fig5_LHR_vs_ONI_range.png"), fig5, width = 7, height = 6, dpi = 300)

cat("\n=== Done. Outputs saved to:", res.dir, "===\n")
