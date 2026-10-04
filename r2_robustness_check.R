###############################################################################
# Signal-strength (R2) robustness check for Arms B and C
# Re-aggregates sim_summary_v8.csv (condition-level) by R2_total; no re-simulation.
# Author: Joseph L. Hagan, ScD, MSPH   Date: October 2026
###############################################################################
library(dplyr)
library(tidyr)

# ── PATHS ────────────────────────────────────────────────────────────────────
# Run from the repository root (setwd() to the cloned mixed-predictors-cv folder),
# where sim_summary_v8.csv is located. Output CSVs are written to results/supplementary/.
IN_FILE <- "sim_summary_v8.csv"
OUT_DIR <- file.path("results", "supplementary")

stopifnot(file.exists(IN_FILE))   # fails with a clear message if run from the wrong folder
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── READ AND SPLIT BY ARM ────────────────────────────────────────────────────
d <- read.csv(IN_FILE, stringsAsFactors = FALSE) %>%
  mutate(Learner = recode(learner, ridge_linear = "LR", xgboost = "XGB"))

dA <- filter(d, arm == "A_fixed")
dB <- filter(d, arm == "B_long")
dC <- filter(d, arm == "C_mixed")

# 1. Arm A: R2 gradient (Section 3.1 sentence)
tabA <- dA %>% group_by(Learner, R2_total) %>%
  summarise(delta_naive = mean(mean_delta_naive), .groups = "drop")

# 2. Arm B: ICC gradient by R2
tabB_icc <- dB %>% group_by(Learner, R2_total, ICC) %>%
  summarise(delta_naive = mean(mean_delta_naive), .groups = "drop")

# 3. Arm B: effect of rho (0.7 vs 0.3), % change, by R2 and ICC
tabB_rho <- dB %>% group_by(Learner, R2_total, ICC, rho) %>%
  summarise(delta = mean(mean_delta_naive), .groups = "drop") %>%
  pivot_wider(names_from = rho, values_from = delta, names_prefix = "rho_") %>%
  mutate(pct_change = 100 * (rho_0.7 / rho_0.3 - 1))

# 4. Arm B (ICC = 0.9) vs Arm A at matched predictor count, n, and R2
a <- dA %>% select(Learner, R2_total, n, p = p_F, delta_A = mean_delta_naive)
b <- dB %>% filter(ICC == 0.9) %>%
  group_by(Learner, R2_total, n, p = p_L) %>%
  summarise(delta_B = mean(mean_delta_naive), .groups = "drop")
tabConv <- inner_join(a, b, by = c("Learner", "R2_total", "n", "p")) %>%
  mutate(diff = delta_B - delta_A)
tabConv_sum <- tabConv %>% group_by(Learner, R2_total) %>%
  summarise(mean_A = mean(delta_A), mean_B = mean(delta_B),
            max_abs_diff = max(abs(diff)), .groups = "drop")

# 5. Arm C: signal allocation by R2
tabC <- dC %>% group_by(Learner, R2_total, signal_allocation) %>%
  summarise(delta_naive = mean(mean_delta_naive), .groups = "drop") %>%
  pivot_wider(names_from = signal_allocation, values_from = delta_naive) %>%
  mutate(long_minus_fixed = long_dominant - fixed_dominant)

# 6. Ceiling check: proportion of conditions with mean naive AUROC >= 0.999
ceil <- d %>% group_by(arm, Learner, R2_total) %>%
  summarise(prop_ceiling = mean(mean_auroc_naive >= 0.999), .groups = "drop")

# ── PRINT AND SAVE ───────────────────────────────────────────────────────────
res <- list(armA_R2 = tabA, armB_ICC_by_R2 = tabB_icc, armB_rho_by_R2 = tabB_rho,
            armB_vs_A_convergence = tabConv_sum, armC_allocation_by_R2 = tabC,
            ceiling = ceil)
for (nm in names(res)) {
  cat("\n", nm, "\n", sep = "")
  print(as.data.frame(res[[nm]]), digits = 3)
  write.csv(res[[nm]], file.path(OUT_DIR, paste0("r2_check_", nm, ".csv")), row.names = FALSE)
}
cat("\nCSV files written to:", OUT_DIR, "\n")
