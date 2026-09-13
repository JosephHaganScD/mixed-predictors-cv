###############################################################################
# supplementary_tables_S1_S10.R
#
# Produces all ten supplementary tables (S1 through S10) for:
#   "Identity-Mediated Leakage in Naive Cross-Validation of Clinical
#    Prediction Models with Fixed and Repeated-Measures Predictors"
#   Target journal: Statistics in Medicine
#
# Output: a single Word document (Supplementary_Tables_S1_S10.docx) written
#   to OUT_DIR, with each table preceded by its title/legend paragraph and
#   separated by a page break.
#
# Data sources (all paths set in Section 1):
#   sim_results_v8.csv    -- iteration-level (567,000 rows); Brier columns
#   sim_summary_v8.csv    -- condition-level summary (1,134 rows)
#   fingerprint_stability_summary.csv
#   mixed_fingerprint_stability.csv
#   categorical_results.csv
#   penalty_results.csv
#
# S4 and S10 require sim_results_v8.csv. If that file is absent the script
# will skip those two tables and warn; all others will still be produced.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Section of Neonatology, Baylor College of Medicine
###############################################################################

library(dplyr)
library(tidyr)
library(officer)
library(flextable)

# =============================================================================
# 1. PATHS
# =============================================================================
# Set these three directories before running:
#
#   SIM_DIR  : folder containing sim_results_v8.csv (iteration-level, ~208 MB).
#              This file is NOT included in the repository because it exceeds
#              GitHub file size limits. Regenerate it with run_simulation.R,
#              or set SIM_DIR <- NULL to skip S4 and S10.
#
#   DATA_DIR : folder containing the five summary-level CSV files
#              (sim_summary_v8.csv, categorical_results.csv, etc.).
#              If you cloned this repository, set DATA_DIR <- "." to use the
#              current working directory.
#
#   OUT_DIR  : folder where Supplementary_Tables_S1_S10.docx will be saved.
#              "." saves to the current working directory.

SIM_DIR  <- "."   # change to folder containing sim_results_v8.csv, or NULL to skip S4/S10
DATA_DIR <- "."   # change to folder containing the summary CSVs (or leave as "." if running from repo)
OUT_DIR  <- "."   # change to desired output folder

SIM_SUMMARY   <- file.path(SIM_DIR,  "sim_summary_v8.csv")
SIM_ITER      <- file.path(SIM_DIR,  "sim_results_v8.csv")   # large; S4 + S10 only
FPRINT_STAB   <- file.path(DATA_DIR, "fingerprint_stability_summary.csv")
MIX_FPRINT    <- file.path(DATA_DIR, "mixed_fingerprint_stability.csv")
CAT_RESULTS   <- file.path(DATA_DIR, "categorical_results.csv")
PENALTY       <- file.path(DATA_DIR, "penalty_results.csv")

# =============================================================================
# 2. HELPERS
# =============================================================================

bp <- function(color = "black", width = 0.5) fp_border(color = color, width = width)

# Standard flextable style (Times New Roman, minimalist borders, 9pt body)
style_ft <- function(ft, hline_rows = NULL) {
  ft <- ft %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 9,  part = "body") %>%
    fontsize(size = 9,  part = "header") %>%
    bold(part = "header") %>%
    align(part = "header", align = "center") %>%
    border_remove() %>%
    hline_top(border = bp(width = 1.0), part = "header") %>%
    hline_bottom(border = bp(width = 1.0), part = "header") %>%
    hline_bottom(border = bp(width = 1.0), part = "body")

  if (!is.null(hline_rows)) {
    for (r in hline_rows)
      ft <- hline(ft, i = r, border = bp(width = 0.5))
  }

  ft %>% autofit() %>% fit_to_width(max_width = 9)
}

# Add a table to the document with title paragraph + optional page break after
add_table <- function(doc, title_text, ft, page_break = TRUE) {
  doc <- body_add_par(doc, title_text, style = "Normal")
  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_flextable(doc, ft)
  if (page_break) doc <- body_add_break(doc)
  doc
}

fmt3  <- function(x) sprintf("%.3f", as.numeric(x))
fmtp  <- function(x) sprintf("%+.3f", as.numeric(x))
fmt4  <- function(x) sprintf("%.4f", as.numeric(x))
fmt5  <- function(x) sprintf("%.5f", as.numeric(x))

# Learner display labels
lrn_label <- function(x) ifelse(x %in% c("ridge_linear", "logistic"), "LR", "XGBoost")

# =============================================================================
# 3. READ DATA
# =============================================================================

cat("Reading sim_summary_v8.csv ...\n")
smry <- read.csv(SIM_SUMMARY, stringsAsFactors = FALSE)

cat("Reading fingerprint_stability_summary.csv ...\n")
fstab <- read.csv(FPRINT_STAB, stringsAsFactors = FALSE)

cat("Reading mixed_fingerprint_stability.csv ...\n")
mstab <- read.csv(MIX_FPRINT, stringsAsFactors = FALSE)

cat("Reading categorical_results.csv ...\n")
cat_res <- read.csv(CAT_RESULTS, stringsAsFactors = FALSE)

cat("Reading penalty_results.csv ...\n")
pen <- read.csv(PENALTY, stringsAsFactors = FALSE)

# Attempt to read iteration-level file (large; needed for S4 and S10)
have_iter <- !is.null(SIM_DIR) && file.exists(SIM_ITER)
if (have_iter) {
  cat("Reading sim_results_v8.csv (this may take a moment) ...\n")
  iter <- read.csv(SIM_ITER, stringsAsFactors = FALSE)
  cat(sprintf("  Read %d rows.\n", nrow(iter)))
} else {
  warning("sim_results_v8.csv not found at ", SIM_ITER,
          "\nS4 and S10 will be skipped.")
}

# =============================================================================
# 4. OPEN DOCUMENT
# =============================================================================

doc <- read_docx()

# =============================================================================
# S1. Full Arm A condition-level table, 27 conditions x 2 learners,
#     with Monte Carlo standard errors
# =============================================================================
cat("Building S1 ...\n")

s1 <- smry %>%
  filter(arm == "A_fixed") %>%
  mutate(
    Learner    = lrn_label(learner),
    n          = as.integer(n),
    `p_F`      = as.integer(p_F),
    R2         = as.numeric(R2_total),
    `AUROC_ext`= fmt3(mean_auroc_true),
    `AUROC naive` = fmt3(mean_delta_naive + mean_auroc_true),
    `AUROC subject` = fmt3(mean_delta_subject + mean_auroc_true),
    `Delta naive`   = fmtp(mean_delta_naive),
    `Delta subject` = fmtp(mean_delta_subject),
    `MCSE naive`    = fmt5(mcse_delta_naive),
    `MCSE subject`  = fmt5(mcse_delta_subject)
  ) %>%
  arrange(Learner, n, p_F, R2) %>%
  select(Learner, n, `p_F`, R2,
         `AUROC_ext`, `AUROC naive`, `AUROC subject`,
         `Delta naive`, `Delta subject`,
         `MCSE naive`, `MCSE subject`)

# Hline after each learner block (27 conditions per learner)
ft_s1 <- flextable(s1) %>%
  merge_v(j = c("Learner", "n", "p_F")) %>%
  valign(j = c("Learner","n","p_F"), valign = "top") %>%
  align(j = 3:11, align = "center", part = "body") %>%
  align(j = 1:2,  align = "left",   part = "body") %>%
  style_ft(hline_rows = 27)

doc <- add_table(doc,
  paste0("Supplementary Table S1. Full Arm A condition-level results, 27 conditions ",
         "(p_F in {2, 5, 10}, n in {50, 100, 200}, R\u00b2 in {0.05, 0.15, 0.30}), ",
         "both learners, 500 iterations per condition. AUROC_ext is the mean ",
         "out-of-sample AUROC on the independent 5,000-subject reference cohort. ",
         "Delta values are optimism relative to AUROC_ext. MCSE, Monte Carlo ",
         "standard error. LR, penalized logistic regression."),
  ft_s1)

# =============================================================================
# S2. Full Arm B condition-level table, collapsed over signal strength (R2)
# =============================================================================
cat("Building S2 ...\n")

s2 <- smry %>%
  filter(arm == "B_long") %>%
  group_by(learner, p_L, ICC, rho, n) %>%
  summarise(
    mean_delta_naive   = mean(mean_delta_naive),
    mean_delta_subject = mean(mean_delta_subject),
    mcse_naive         = mean(mcse_delta_naive),
    mcse_subject       = mean(mcse_delta_subject),
    .groups = "drop"
  ) %>%
  mutate(
    Learner  = lrn_label(learner),
    p_L      = as.integer(p_L),
    ICC      = as.numeric(ICC),
    rho      = as.numeric(rho),
    n        = as.integer(n),
    `Delta naive`    = fmtp(mean_delta_naive),
    `Delta subject`  = fmtp(mean_delta_subject),
    `MCSE naive`     = fmt5(mcse_naive),
    `MCSE subject`   = fmt5(mcse_subject)
  ) %>%
  arrange(Learner, p_L, ICC, rho, n) %>%
  select(Learner, p_L, ICC, rho, n,
         `Delta naive`, `Delta subject`,
         `MCSE naive`, `MCSE subject`)

# Hline between learner blocks
n_b_per_learner <- nrow(s2) / 2

ft_s2 <- flextable(s2) %>%
  merge_v(j = c("Learner","p_L","ICC","rho")) %>%
  valign(j = c("Learner","p_L","ICC","rho"), valign = "top") %>%
  align(j = 3:9, align = "center", part = "body") %>%
  align(j = 1:2, align = "left",   part = "body") %>%
  style_ft(hline_rows = n_b_per_learner)

doc <- add_table(doc,
  paste0("Supplementary Table S2. Full Arm B condition-level results, 108 conditions ",
         "(p_L in {2, 5}, ICC in {0.3, 0.7, 0.9}, \u03c1 in {0.3, 0.7}, ",
         "n in {50, 100, 200}), both learners, averaged over signal strength ",
         "(R\u00b2 in {0.05, 0.15, 0.30}), 500 iterations per condition. ",
         "Delta values are mean naive and subject-level optimism relative to ",
         "AUROC_ext. MCSE, Monte Carlo standard error (mean across R\u00b2 levels). ",
         "LR, penalized logistic regression."),
  ft_s2)

# =============================================================================
# S3. Full Arm C condition-level table, collapsed over dependence parameters
#     (ICC, rho) and signal strength (R2)
# =============================================================================
cat("Building S3 ...\n")

s3 <- smry %>%
  filter(arm == "C_mixed") %>%
  group_by(learner, p_F, p_L, n) %>%
  summarise(
    mean_delta_naive   = mean(mean_delta_naive),
    mean_delta_subject = mean(mean_delta_subject),
    mcse_naive         = mean(mcse_delta_naive),
    mcse_subject       = mean(mcse_delta_subject),
    .groups = "drop"
  ) %>%
  mutate(
    Learner = lrn_label(learner),
    p_F     = as.integer(p_F),
    p_L     = as.integer(p_L),
    n       = as.integer(n),
    `Delta naive`   = fmtp(mean_delta_naive),
    `Delta subject` = fmtp(mean_delta_subject),
    `MCSE naive`    = fmt5(mcse_naive),
    `MCSE subject`  = fmt5(mcse_subject)
  ) %>%
  arrange(Learner, p_F, p_L, n) %>%
  select(Learner, p_F, p_L, n,
         `Delta naive`, `Delta subject`,
         `MCSE naive`, `MCSE subject`)

n_c_per_learner <- nrow(s3) / 2

ft_s3 <- flextable(s3) %>%
  merge_v(j = c("Learner","p_F","p_L")) %>%
  valign(j = c("Learner","p_F","p_L"), valign = "top") %>%
  align(j = 3:8, align = "center", part = "body") %>%
  align(j = 1:2, align = "left",   part = "body") %>%
  style_ft(hline_rows = n_c_per_learner)

doc <- add_table(doc,
  paste0("Supplementary Table S3. Full Arm C condition-level results, 144 conditions ",
         "(p_F in {2, 5, 10}, p_L in {2, 5}, n in {50, 100, 200}), both learners, ",
         "averaged over dependence parameters (ICC in {0.3, 0.7}, \u03c1 in {0.3, 0.7}) ",
         "and signal strength (R\u00b2 in {0.05, 0.15, 0.30}), 500 iterations per condition. ",
         "Delta values are mean naive and subject-level optimism relative to AUROC_ext. ",
         "MCSE, Monte Carlo standard error (mean across collapsed levels). ",
         "LR, penalized logistic regression."),
  ft_s3)

# =============================================================================
# S4. Brier-scale optimism for selected contrasts (Arm A predictor-count
#     gradient). Requires sim_results_v8.csv.
# =============================================================================
cat("Building S4 ...\n")

if (have_iter) {

  # Compute per-condition mean Brier optimism from iteration-level data
  # brier_true = Brier_ext (full-model on independent cohort)
  # brier_naive = naive CV Brier
  # brier_subject = subject-level CV Brier
  # Optimism = brier_naive - brier_true (positive = optimism, consistent
  #   with AUROC convention where higher CV AUROC > true AUROC)
  # Note: for Brier, lower is better, so naive < true means the CV
  #   underestimates error, i.e., appears better than it is.
  #   Optimism = brier_true - brier_naive (positive means CV too optimistic)

  brier_smry <- iter %>%
    group_by(condition_id, arm, learner, p_F, p_L, n, R2 = R2_total,
             ICC, rho) %>%
    summarise(
      mean_brier_opt_naive   = mean(brier_true - brier_naive,   na.rm = TRUE),
      mean_brier_opt_subject = mean(brier_true - brier_subject, na.rm = TRUE),
      mcse_brier_naive       = sd(brier_true - brier_naive,     na.rm = TRUE) /
                                 sqrt(n()),
      .groups = "drop"
    )

  # S4: Arm A only, showing predictor-count gradient by n and R2
  s4 <- brier_smry %>%
    filter(arm == "A_fixed") %>%
    mutate(
      Learner = lrn_label(learner),
      p_F     = as.integer(p_F),
      n       = as.integer(n),
      R2      = as.numeric(R2),
      `Brier Delta naive`   = fmtp(mean_brier_opt_naive),
      `Brier Delta subject` = fmtp(mean_brier_opt_subject),
      `MCSE`                = fmt5(mcse_brier_naive)
    ) %>%
    arrange(Learner, n, p_F, R2) %>%
    select(Learner, n, p_F, R2,
           `Brier Delta naive`, `Brier Delta subject`, `MCSE`)

  ft_s4 <- flextable(s4) %>%
    merge_v(j = c("Learner","n","p_F")) %>%
    valign(j = c("Learner","n","p_F"), valign = "top") %>%
    align(j = 3:7, align = "center", part = "body") %>%
    align(j = 1:2, align = "left",   part = "body") %>%
    style_ft(hline_rows = 27)

  doc <- add_table(doc,
    paste0("Supplementary Table S4. Brier-scale optimism for Arm A (fixed-only ",
           "predictors), 27 conditions, both learners. Brier Delta naive is the ",
           "mean of (Brier_ext \u2212 Brier_naive) across 500 iterations; positive ",
           "values indicate optimism. The predictor-count gradient is recovered ",
           "on this scale in conditions where AUROC naive is censored near 1.000. ",
           "MCSE, Monte Carlo standard error of the naive Brier optimism. ",
           "LR, penalized logistic regression."),
    ft_s4)

} else {
  doc <- body_add_par(doc,
    "Supplementary Table S4. [SKIPPED: sim_results_v8.csv not found]",
    style = "Normal")
  doc <- body_add_break(doc)
}

# =============================================================================
# S5. Fingerprint stability across the full ICC grid, with uniqueness shown
# =============================================================================
cat("Building S5 ...\n")

s5 <- fstab %>%
  mutate(
    ICC       = as.numeric(ICC),
    rho       = as.numeric(rho),
    p_L       = as.integer(p_L),
    n         = as.integer(n),
    `u (Q4)`  = fmt3(u_q4),
    `u (Q10)` = fmt3(u_q10),
    `Stability (Q4)`  = fmt3(stab_q4),
    `Stability (Q10)` = fmt3(stab_q10),
    `MCSE stab (Q4)`  = fmt5(mcse_stab_q4),
    `MCSE stab (Q10)` = fmt5(mcse_stab_q10)
  ) %>%
  arrange(p_L, rho, n, ICC) %>%
  select(p_L, rho, n, ICC,
         `u (Q4)`, `u (Q10)`,
         `Stability (Q4)`, `Stability (Q10)`,
         `MCSE stab (Q4)`, `MCSE stab (Q10)`)

ft_s5 <- flextable(s5) %>%
  merge_v(j = c("p_L","rho","n")) %>%
  valign(j = c("p_L","rho","n"), valign = "top") %>%
  align(j = 3:10, align = "center", part = "body") %>%
  align(j = 1:2,  align = "left",   part = "body") %>%
  style_ft()

doc <- add_table(doc,
  paste0("Supplementary Table S5. Fingerprint stability and subject uniqueness ",
         "across the full intraclass correlation grid, 108 conditions ",
         "(p_L in {2, 5}, ICC in {0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1}, ",
         "\u03c1 in {0.3, 0.7}, n in {50, 100, 200}). u (Q4) and u (Q10) are the mean ",
         "proportion of subjects occupying a singleton equivalence class under ",
         "quartile and decile binning, respectively. Stability is the mean ",
         "fingerprint stability statistic. MCSE, Monte Carlo standard error."),
  ft_s5)

# =============================================================================
# S6. Composition factorial, full breakdown by sample size and signal strength
# =============================================================================
cat("Building S6 ...\n")

# Ceiling: AUROC naive >= 0.999
s6 <- cat_res %>%
  mutate(
    Learner     = lrn_label(learner),
    Composition = composition,
    n           = as.integer(n),
    R2          = as.numeric(R2),
    u_raw_fmt   = sprintf("%.3f", as.numeric(u_raw)),
    auroc_naive_fmt   = fmt3(auroc_naive.auroc),
    auroc_subject_fmt = fmt3(auroc_subject.auroc),
    delta_naive_fmt   = fmtp(as.numeric(auroc_naive.auroc) -
                               as.numeric(auroc_true.auroc)),
    ceiling_pct = ifelse(as.numeric(auroc_naive.auroc) >= 0.999, "Yes", "No")
  ) %>%
  arrange(Learner,
          factor(Composition, levels = c("5 continuous","3 continuous + 2 binary",
                                         "2 continuous + 3 binary",
                                         "1 continuous + 4 binary","5 binary")),
          n, R2) %>%
  select(Learner, Composition, n, R2, u_raw_fmt,
         auroc_naive_fmt, auroc_subject_fmt, delta_naive_fmt, ceiling_pct)

names(s6) <- c("Learner","Composition","n","R\u00b2","u",
               "AUROC naive","AUROC subject","Delta naive",
               "At ceiling")

n_cat_per_learner <- nrow(s6) / 2

ft_s6 <- flextable(s6) %>%
  merge_v(j = c("Learner","Composition","n")) %>%
  valign(j = c("Learner","Composition","n"), valign = "top") %>%
  align(j = 3:9, align = "center", part = "body") %>%
  align(j = 1:2, align = "left",   part = "body") %>%
  style_ft(hline_rows = n_cat_per_learner)

doc <- add_table(doc,
  paste0("Supplementary Table S6. Composition factorial, full breakdown by sample ",
         "size and signal strength, all five predictor compositions, both learners. ",
         "u is the mean raw subject uniqueness proportion. AUROC naive and AUROC ",
         "subject are means across 500 iterations. Delta naive is optimism relative ",
         "to AUROC_ext. At ceiling indicates whether mean naive AUROC reached or ",
         "exceeded 0.999. LR, penalized logistic regression."),
  ft_s6)

# =============================================================================
# S7. Penalty sensitivity: coefficient norms and optimism at all six penalty
#     levels, 12 conditions
# =============================================================================
cat("Building S7 ...\n")

# Summarise over iterations within condition x lambda
s7_raw <- pen %>%
  mutate(lambda_num = as.numeric(lambda)) %>%
  group_by(condition_id, arm, p_F, p_L, n, R2, lambda_num, is_main_rule) %>%
  summarise(
    mean_coef_l2     = mean(as.numeric(coef_l2),     na.rm = TRUE),
    mean_auroc_naive = mean(as.numeric(auroc_naive),  na.rm = TRUE),
    mean_auroc_subj  = mean(as.numeric(auroc_subject),na.rm = TRUE),
    mean_auroc_true  = mean(as.numeric(auroc_true),   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    delta_naive   = mean_auroc_naive - mean_auroc_true,
    delta_subject = mean_auroc_subj  - mean_auroc_true,
    Arm           = ifelse(arm == "A_fixed", "A (fixed)", "B (long)"),
    p_label       = ifelse(arm == "A_fixed",
                           paste0("p_F=", p_F),
                           paste0("p_L=", p_L)),
    Condition     = paste0(Arm, ", ", p_label, ", n=", n, ", R\u00b2=", R2),
    Lambda        = sprintf("%.2e", lambda_num),
    `Main rule`   = ifelse(is_main_rule == "TRUE", "Yes", "No"),
    `||beta||`    = fmt3(mean_coef_l2),
    `AUROC naive` = fmt3(mean_auroc_naive),
    `Delta naive` = fmtp(delta_naive),
    `Delta subject` = fmtp(delta_subject)
  ) %>%
  arrange(as.integer(condition_id), lambda_num) %>%
  select(Condition, Lambda, `Main rule`, `||beta||`,
         `AUROC naive`, `Delta naive`, `Delta subject`)

# Hline between conditions (6 lambda levels per condition)
hlines_s7 <- seq(6, 6 * 11, by = 6)

ft_s7 <- flextable(s7_raw) %>%
  merge_v(j = "Condition") %>%
  valign(j = "Condition", valign = "top") %>%
  align(j = 2:7, align = "center", part = "body") %>%
  align(j = 1,   align = "left",   part = "body") %>%
  style_ft(hline_rows = hlines_s7)

doc <- add_table(doc,
  paste0("Supplementary Table S7. Penalty sensitivity analysis. Six penalty levels ",
         "evaluated across 12 conditions spanning the extremes of Arms A and B. ",
         "||beta|| is the mean L2 coefficient norm. Main rule identifies the ",
         "penalty used in the primary analysis (\u03bb = log(1+1/n)/n). ",
         "Delta values are optimism relative to AUROC_ext. ",
         "LR only (penalized logistic regression)."),
  ft_s7)

# =============================================================================
# S8. Mixed-fingerprint stability by composition, including complementary
#     disambiguation result for categorical fixed sets
# =============================================================================
cat("Building S8 ...\n")

s8 <- mstab %>%
  mutate(
    p_L   = as.integer(p_L),
    ICC   = as.numeric(ICC),
    rho   = as.numeric(rho),
    n     = as.integer(n),
    Composition = fixed_comp,
    `s_F`     = ifelse(is.na(s_F),    "—", fmt3(s_F)),
    `s_L`     = fmt3(s_L),
    `s_joint` = fmt3(s_joint),
    `Increment` = ifelse(is.na(increment), "—",
                         sprintf("%+.3f", as.numeric(increment)))
  ) %>%
  arrange(factor(Composition,
                 levels = c("none","1 continuous","2 continuous",
                            "2 continuous + 3 binary","5 continuous","5 binary")),
          p_L, ICC, rho, n) %>%
  select(Composition, p_L, ICC, rho, n,
         `s_F`, `s_L`, `s_joint`, `Increment`)

ft_s8 <- flextable(s8) %>%
  merge_v(j = c("Composition","p_L","ICC","rho")) %>%
  valign(j = c("Composition","p_L","ICC","rho"), valign = "top") %>%
  align(j = 4:9, align = "center", part = "body") %>%
  align(j = 1:3, align = "left",   part = "body") %>%
  style_ft()

doc <- add_table(doc,
  paste0("Supplementary Table S8. Mixed-fingerprint stability by fixed-predictor ",
         "composition, 216 conditions. s_F is the fixed-fingerprint stability; ",
         "s_L is the longitudinal-fingerprint stability; s_joint is the joint ",
         "fingerprint stability. Increment is s_joint \u2212 max(s_F, s_L), ",
         "indicating the disambiguation contribution of combining fixed and ",
         "longitudinal fingerprints; a negative increment indicates that the ",
         "joint fingerprint is less stable than the better of the two individual ",
         "fingerprints (complementary disambiguation). \u2014 denotes conditions ",
         "with no fixed predictors. p_L, number of longitudinal predictors; ",
         "ICC, intraclass correlation; \u03c1, lag-1 autocorrelation."),
  ft_s8)

# =============================================================================
# S9. Empirical illustration: secondary outcome and reduced-predictor
#     sensitivity analysis (V2: Avg_FiO2 + Avg_SpO2)
#     All values hard-coded from verified analysis output files.
# =============================================================================
cat("Building S9 ...\n")

# --- Panel A: Secondary outcome (severe ROP) ---
# Source: rop_decomposition_xgb_2026-07-02.txt
panel_a <- data.frame(
  Learner       = c("LR","LR","LR","XGBoost","XGBoost","XGBoost"),
  Configuration = rep(c("Fixed-only","Longitudinal-only","Mixed"), 2),
  auroc_naive   = c(0.864, 0.700, 0.885, 1.000, 0.793, 1.000),
  auroc_subject = c(0.706, 0.532, 0.708, 0.735, 0.562, 0.747),
  auroc_loco    = c(0.701, 0.507, 0.707, 0.760, 0.566, 0.755),
  loco_lo       = c(0.507, 0.365, 0.513, 0.627, 0.431, 0.613),
  loco_hi       = c(0.896, 0.648, 0.902, 0.893, 0.702, 0.896),
  delta_naive   = c(+0.163, +0.193, +0.177, +0.240, +0.227, +0.245),
  delta_subject = c(+0.005, +0.025, +0.001, -0.025, -0.004, -0.008),
  stringsAsFactors = FALSE
)

# --- Panel B: Reduced predictor set, V2, primary outcome ---
# Source: rop_tier2_sensitivity_2026-07-02.txt (V2 block)
panel_b <- data.frame(
  Learner       = c("LR","LR","LR","XGBoost","XGBoost","XGBoost"),
  Configuration = rep(c("Fixed-only","Longitudinal-only","Mixed"), 2),
  auroc_naive   = c(0.657, 0.667, 0.688, 0.999, 0.688, 1.000),
  auroc_subject = c(0.607, 0.587, 0.624, 0.542, 0.603, 0.547),
  auroc_loco    = c(0.609, 0.603, 0.624, 0.557, 0.616, 0.557),
  loco_lo       = c(0.497, 0.491, 0.513, 0.443, 0.504, 0.443),
  loco_hi       = c(0.722, 0.715, 0.735, 0.670, 0.727, 0.672),
  delta_naive   = c(+0.048, +0.064, +0.065, +0.442, +0.073, +0.442),
  delta_subject = c(-0.002, -0.016, +0.001, -0.014, -0.013, -0.010),
  stringsAsFactors = FALSE
)

fmt_panel <- function(df) {
  df %>%
    mutate(
      `AUROC, naive`         = fmt3(auroc_naive),
      `AUROC, subject-level` = fmt3(auroc_subject),
      `AUROC, LOCO (95% CI)` = paste0(fmt3(auroc_loco),
                                       " (", fmt3(loco_lo), " to ",
                                       fmt3(loco_hi), ")"),
      `Delta naive`          = fmtp(delta_naive),
      `Delta subject-level`  = fmtp(delta_subject)
    ) %>%
    select(Learner, Configuration,
           `AUROC, naive`, `AUROC, subject-level`, `AUROC, LOCO (95% CI)`,
           `Delta naive`, `Delta subject-level`)
}

make_s9_ft <- function(df) {
  flextable(df) %>%
    merge_v(j = "Learner") %>%
    valign(j = "Learner", valign = "top") %>%
    align(j = 1:2, align = "left",   part = "body") %>%
    align(j = 3:7, align = "center", part = "body") %>%
    style_ft(hline_rows = 3)
}

ft_s9a <- make_s9_ft(fmt_panel(panel_a))
ft_s9b <- make_s9_ft(fmt_panel(panel_b))

# Add S9 as two sub-tables within one page section
s9_title <- paste0(
  "Supplementary Table S9. Empirical illustration: secondary outcome and ",
  "reduced-predictor sensitivity analysis. AUROC values are means across 50 ",
  "fold-assignment replications for naive and subject-level cross-validation; ",
  "the leave-one-cluster-out (LOCO) estimate is deterministic with a 95% ",
  "confidence interval by the DeLong method. Delta values are optimism relative ",
  "to LOCO. LR, penalized logistic regression; LOCO, leave-one-cluster-out."
)

doc <- body_add_par(doc, s9_title, style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc,
  paste0("Panel A. Secondary outcome (severe retinopathy of prematurity, ",
         "15 events among 101 subjects). Full predictor set: 8 fixed covariates ",
         "and 7 longitudinal oxygenation variables."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_flextable(doc, ft_s9a)
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_par(doc,
  paste0("Panel B. Reduced predictor set matched to simulation conditions ",
         "(p_F = 2, p_L = 2): fixed covariates limited to birth weight and ",
         "gestational age; longitudinal predictors limited to Avg_FiO2 and ",
         "Avg_SpO2. Primary outcome (any retinopathy of prematurity, ",
         "48 events among 101 subjects)."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")
doc <- body_add_flextable(doc, ft_s9b)
doc <- body_add_break(doc)

# =============================================================================
# S10. Condition-level Brier-scale optimism for all three arms.
#      Requires sim_results_v8.csv.
# =============================================================================
cat("Building S10 ...\n")

if (have_iter) {

  if (!exists("brier_smry")) {
    brier_smry <- iter %>%
      group_by(condition_id, arm, learner, p_F, p_L, n,
               R2 = R2_total, ICC, rho) %>%
      summarise(
        mean_brier_opt_naive   = mean(brier_true - brier_naive,   na.rm = TRUE),
        mean_brier_opt_subject = mean(brier_true - brier_subject, na.rm = TRUE),
        mcse_brier_naive       = sd(brier_true - brier_naive,     na.rm = TRUE) /
                                   sqrt(n()),
        .groups = "drop"
      )
  }

  # S10: all three arms, collapsed the same way as S1-S3
  # Arm A: no collapsing needed (27 conditions)
  # Arm B: collapse over R2
  # Arm C: collapse over ICC, rho, R2

  s10_a <- brier_smry %>%
    filter(arm == "A_fixed") %>%
    mutate(Arm = "A", Learner = lrn_label(learner),
           p_F = as.integer(p_F), n = as.integer(n), R2 = as.numeric(R2),
           Condition = paste0("p_F=", p_F, ", n=", n, ", R\u00b2=", R2)) %>%
    arrange(Learner, n, p_F, R2)

  s10_b <- brier_smry %>%
    filter(arm == "B_long") %>%
    group_by(learner, p_L, ICC, rho, n) %>%
    summarise(
      mean_brier_opt_naive   = mean(mean_brier_opt_naive),
      mean_brier_opt_subject = mean(mean_brier_opt_subject),
      mcse_brier_naive       = mean(mcse_brier_naive),
      .groups = "drop"
    ) %>%
    mutate(Arm = "B", Learner = lrn_label(learner),
           p_L = as.integer(p_L), ICC = as.numeric(ICC),
           rho = as.numeric(rho), n = as.integer(n),
           Condition = paste0("p_L=", p_L, ", ICC=", ICC,
                              ", \u03c1=", rho, ", n=", n)) %>%
    arrange(Learner, p_L, ICC, rho, n)

  s10_c <- brier_smry %>%
    filter(arm == "C_mixed") %>%
    group_by(learner, p_F, p_L, n) %>%
    summarise(
      mean_brier_opt_naive   = mean(mean_brier_opt_naive),
      mean_brier_opt_subject = mean(mean_brier_opt_subject),
      mcse_brier_naive       = mean(mcse_brier_naive),
      .groups = "drop"
    ) %>%
    mutate(Arm = "C", Learner = lrn_label(learner),
           p_F = as.integer(p_F), p_L = as.integer(p_L),
           n = as.integer(n),
           Condition = paste0("p_F=", p_F, ", p_L=", p_L, ", n=", n)) %>%
    arrange(Learner, p_F, p_L, n)

  # Bind all three arms
  s10 <- bind_rows(
    s10_a %>% select(Arm, Learner, Condition,
                     mean_brier_opt_naive, mean_brier_opt_subject,
                     mcse_brier_naive),
    s10_b %>% select(Arm, Learner, Condition,
                     mean_brier_opt_naive, mean_brier_opt_subject,
                     mcse_brier_naive),
    s10_c %>% select(Arm, Learner, Condition,
                     mean_brier_opt_naive, mean_brier_opt_subject,
                     mcse_brier_naive)
  ) %>%
    mutate(
      `Brier Delta naive`   = fmtp(mean_brier_opt_naive),
      `Brier Delta subject` = fmtp(mean_brier_opt_subject),
      `MCSE`                = fmt5(mcse_brier_naive)
    ) %>%
    select(Arm, Learner, Condition,
           `Brier Delta naive`, `Brier Delta subject`, `MCSE`)

  # Hlines between arms (within each learner block)
  n_a <- nrow(s10_a) / 2
  n_b <- nrow(s10_b) / 2
  n_c <- nrow(s10_c) / 2
  total_per_learner <- n_a + n_b + n_c
  hlines_s10 <- c(n_a, n_a + n_b,                        # within LR block
                  total_per_learner,                       # between learners
                  total_per_learner + n_a,
                  total_per_learner + n_a + n_b)

  ft_s10 <- flextable(s10) %>%
    merge_v(j = c("Arm","Learner")) %>%
    valign(j = c("Arm","Learner"), valign = "top") %>%
    align(j = 3:6, align = "center", part = "body") %>%
    align(j = 1:2, align = "left",   part = "body") %>%
    style_ft(hline_rows = hlines_s10)

  doc <- add_table(doc,
    paste0("Supplementary Table S10. Condition-level Brier-scale optimism for ",
           "all three simulation arms, both learners. Brier Delta naive is the ",
           "mean of (Brier_ext \u2212 Brier_naive) across 500 iterations per condition; ",
           "Brier Delta subject is the corresponding subject-level quantity. ",
           "Positive values indicate optimism (naive CV underestimates prediction ",
           "error). Arms B and C are collapsed over signal strength; Arm C is ",
           "additionally collapsed over dependence parameters, consistent with ",
           "Tables S2 and S3. MCSE, Monte Carlo standard error. ",
           "LR, penalized logistic regression."),
    ft_s10, page_break = FALSE)

} else {
  doc <- body_add_par(doc,
    "Supplementary Table S10. [SKIPPED: sim_results_v8.csv not found]",
    style = "Normal")
}

# =============================================================================
# 5. SAVE
# =============================================================================

out_path <- file.path(OUT_DIR, "Supplementary_Tables_S1_S10.docx")
print(doc, target = out_path)
cat(sprintf("\nDone. Output written to:\n  %s\n", out_path))
