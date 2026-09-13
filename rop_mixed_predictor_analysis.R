# Set the path to your copy of the ROP dataset (final daily data.csv).
# Note: this dataset is not publicly available (see README for details).
rop <- read.csv("final daily data.csv")   # update path as needed
head(rop)
names(rop)

###############################################################################
# rop_mixed_predictor_analysis.R
# Empirical Illustration: Mixed-Predictor CV Optimism Decomposition
# ROP Dataset (Srivatsa, Northside Hospital, 2016-2020)
#
# Purpose: Demonstrate the leakage decomposition framework using the ROP
# dataset as a motivating example. Three predictor configurations are
# evaluated, corresponding to the three arms of the simulation study:
#
#   Arm A (fixed-only):      BWGT, GA, AP5, black, white, asian, other, HISP
#   Arm B (longitudinal):    7 oxygenation variables (daily, long format)
#   Arm C (mixed):           Fixed + longitudinal
#
# For each configuration, three CV strategies are applied:
#   1. Naive 10-fold CV (observation-level fold assignment)
#   2. Subject-level 10-fold CV (cluster-aware)
#   3. LOCO CV (leave-one-cluster-out; feasible at n=101)
#
# Optimism = AUROC_strategy - AUROC_LOCO (LOCO serves as the reference)
# Decomposition test: Delta_naive,mixed ≈ Delta_naive,fixed + Delta_naive,long
#
# Primary outcome:   any_rop  (ROP >= 1; 48/101 = 47.5% prevalence)
# Secondary outcome: severe_rop (ROP == 2; 15/101 = 14.9% prevalence)
#
# Learner: Ridge logistic regression with fixed lambda = log(1 + 1/n) / n,
# consistent with the simulation study (Hagan, manuscript under review).
# All models fit in long format for longitudinal and mixed configurations;
# subject-level predicted probabilities obtained by averaging row-level
# predictions within subject.
#
# Author: Joseph L. Hagan, ScD, MSPH
# Section of Neonatology, Baylor College of Medicine
###############################################################################

library(glmnet)
library(pROC)

# =============================================================================
# USER SETTINGS
# =============================================================================

DATA_FILE  <- "data/final_daily_data.csv" # kept for reference only
raw <- rop                                   # use already-loaded object


# Set OUTPUT_DIR to the folder where results should be saved.
OUTPUT_DIR <- "."   # change to your preferred output directory


N_REPS     <- 50    # CV fold assignment replications (for naive and subject CV)
K_FOLDS    <- 10    # number of CV folds

if (!file.exists(DATA_FILE))
  stop(sprintf("Data file not found: %s", DATA_FILE))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

output_file <- file.path(OUTPUT_DIR,
  paste0("rop_decomposition_", format(Sys.Date(), "%Y-%m-%d"), ".txt"))
sink(output_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat(sprintf("ROP Mixed-Predictor Decomposition Analysis\n%s\n\n", Sys.time()))

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

raw <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

# Longitudinal predictor set (7 oxygenation variables, daily)
long_vars <- c("Avg_FiO2", "Avg_SpO2", "Severe_Hypoxemia__",
               "Amb_Hyperox_", "Iatr_HyperOx_", "Swings", "Tit_Index")

# Fixed predictor set (time-invariant subject characteristics)
# Note: sex is not available in this dataset; see manuscript Limitations.
fixed_vars <- c("BWGT", "GA", "AP5", "HISP", "black", "white", "asian", "other")

all_vars <- c(fixed_vars, long_vars)

# Complete-case restriction
cc  <- complete.cases(raw[, all_vars])
dat <- raw[cc, ]
cat(sprintf("Rows after complete-case restriction: %d (%d excluded)\n",
            nrow(dat), sum(!cc)))
cat(sprintf("Subjects: %d\n", length(unique(dat$Case_ID))))
cat(sprintf("Median observations per subject: %g [range %d-%d]\n",
            median(table(dat$Case_ID)),
            min(table(dat$Case_ID)),
            max(table(dat$Case_ID))))

# Binary outcomes
dat$any_rop    <- as.integer(dat$ROP >= 1)
dat$severe_rop <- as.integer(dat$ROP == 2)

n_subj <- length(unique(dat$Case_ID))
cat(sprintf("\nPrimary outcome (any ROP):    %d events / %d subjects (%.1f%%)\n",
            sum(tapply(dat$any_rop,    dat$Case_ID, max)),
            n_subj,
            100 * mean(tapply(dat$any_rop, dat$Case_ID, max))))
cat(sprintf("Secondary outcome (severe ROP): %d events / %d subjects (%.1f%%)\n",
            sum(tapply(dat$severe_rop, dat$Case_ID, max)),
            n_subj,
            100 * mean(tapply(dat$severe_rop, dat$Case_ID, max))))

# =============================================================================
# 2. UNIQUENESS STATISTICS (fixed predictors as quasi-identifiers)
# =============================================================================

cat("\n========== SUBJECT UNIQUENESS (FIXED PREDICTORS) ==========\n")

subj_data <- dat[!duplicated(dat$Case_ID), fixed_vars, drop = FALSE]

# Discretize continuous fixed predictors into quartile bins
bin_col <- function(x, n_bins) {
  cuts <- quantile(x, probs = seq(0, 1, length.out = n_bins + 1), type = 7)
  cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
  as.integer(cut(x, breaks = cuts, labels = FALSE, include.lowest = TRUE))
}

# Quartile binning (Q4) and decile binning (Q10)
for (nb in c(4, 10)) {
  binned <- as.data.frame(lapply(subj_data, function(col) {
    if (length(unique(col)) > nb) bin_col(col, nb) else col
  }))
  key        <- do.call(paste, c(binned, sep = "_"))
  class_sz   <- table(key)
  pct_unique <- mean(class_sz == 1L)
  cat(sprintf("Q%d binning: %.1f%% of subjects uniquely identified",
              nb, 100 * pct_unique))
  cat(sprintf(" (mean class size %.2f, max %d)\n",
              mean(class_sz), max(class_sz)))
}

# Raw continuous: proportion with unique BWGT x GA combination
key_raw    <- paste(subj_data$BWGT, subj_data$GA, sep = "_")
cat(sprintf("Raw BWGT x GA: %.1f%% of subjects uniquely identified\n",
            100 * mean(table(key_raw) == 1L)))

# =============================================================================
# 3. DEPENDENCE STRUCTURE (longitudinal predictors)
# =============================================================================

cat("\n========== DEPENDENCE STRUCTURE (LONGITUDINAL PREDICTORS) ==========\n")

n_bar <- mean(table(dat$Case_ID))

compute_icc <- function(x, grp) {
  fit  <- aov(x ~ factor(grp))
  ms   <- summary(fit)[[1L]][, "Mean Sq"]
  ms_b <- ms[1L]; ms_w <- ms[2L]
  (ms_b - ms_w) / (ms_b + (n_bar - 1) * ms_w)
}

compute_lag1 <- function(v) {
  r <- sapply(unique(dat$Case_ID), function(id) {
    x <- dat[[v]][dat$Case_ID == id]
    if (length(x) < 3L) return(NA_real_)
    cor(x[-length(x)], x[-1L], use = "complete.obs")
  })
  c(median = median(r, na.rm = TRUE),
    q25    = unname(quantile(r, 0.25, na.rm = TRUE)),
    q75    = unname(quantile(r, 0.75, na.rm = TRUE)))
}

cat(sprintf("\n%-22s  %7s  %7s  %6s  %s\n",
            "Variable", "Mean", "SD", "ICC", "Lag-1 rho [IQR]"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (v in long_vars) {
  icc  <- compute_icc(dat[[v]], dat$Case_ID)
  lag1 <- compute_lag1(v)
  cat(sprintf("%-22s  %7.3f  %7.3f  %6.2f  %.2f [%.2f, %.2f]\n",
              v, mean(dat[[v]], na.rm = TRUE), sd(dat[[v]], na.rm = TRUE),
              icc, lag1["median"], lag1["q25"], lag1["q75"]))
}

# =============================================================================
# 4. MODEL AND CV FUNCTIONS
# =============================================================================

# Fixed lambda consistent with simulation: lambda = log(1 + 1/n) / n
# n = number of subjects (not rows)
fixed_lambda <- function(n_subj) log(1 + 1 / n_subj) / n_subj

# Ridge logistic regression: fit on training data, predict on test data
# Operates on long-format matrices; subject-level predictions computed downstream
fit_ridge <- function(X_train, y_train, X_test, lam) {
  fit <- tryCatch(
    glmnet(X_train, y_train, family = "binomial", alpha = 0,
           lambda = lam, standardize = TRUE),
    error = function(e) {
      glmnet(X_train, y_train, family = "binomial", alpha = 0,
             lambda = max(lam, 0.01), standardize = TRUE)
    }
  )
  as.numeric(predict(fit, newx = X_test, s = fit$lambda[1L],
                     type = "response"))
}

# Aggregate row-level predictions to subject level (mean within subject)
# then compute AUROC against subject-level outcome
row_to_subj_auroc <- function(case_ids, row_preds, row_y) {
  subj_p <- tapply(row_preds, case_ids, mean)
  subj_y <- tapply(row_y,    case_ids, function(x) x[1L])
  subj_y <- subj_y[names(subj_p)]
  if (length(unique(subj_y)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(as.numeric(subj_y),
                                   as.numeric(subj_p), quiet = TRUE)))
}

row_to_subj_roc <- function(case_ids, row_preds, row_y) {
  subj_p <- tapply(row_preds, case_ids, mean)
  subj_y <- tapply(row_y,    case_ids, function(x) x[1L])
  subj_y <- subj_y[names(subj_p)]
  if (length(unique(subj_y)) < 2L) return(NULL)
  pROC::roc(as.numeric(subj_y), as.numeric(subj_p), quiet = TRUE)
}

# ---------------------------------------------------------------------------
# Naive 10-fold CV (observation-level fold assignment)
# ---------------------------------------------------------------------------
cv_naive <- function(X, y, case_ids, k = K_FOLDS, seed = 1L) {
  set.seed(seed)
  n       <- nrow(X)
  lam     <- fixed_lambda(length(unique(case_ids)))
  fold_id <- sample(rep(seq_len(k), length.out = n))
  preds   <- rep(NA_real_, n)

  for (f in seq_len(k)) {
    tr <- fold_id != f
    te <- fold_id == f
    if (length(unique(y[tr])) < 2L) next
    preds[te] <- fit_ridge(X[tr, , drop = FALSE], y[tr],
                           X[te, , drop = FALSE], lam)
  }
  valid <- !is.na(preds)
  row_to_subj_auroc(case_ids[valid], preds[valid], y[valid])
}

# ---------------------------------------------------------------------------
# Subject-level 10-fold CV (cluster-aware)
# ---------------------------------------------------------------------------
cv_subject <- function(X, y, case_ids, k = K_FOLDS, seed = 1L) {
  set.seed(seed)
  lam      <- fixed_lambda(length(unique(case_ids)))
  subj_ids <- unique(case_ids)
  n_s      <- length(subj_ids)
  k_act    <- min(k, n_s)
  fold_map <- setNames(rep(seq_len(k_act), length.out = n_s),
                       subj_ids[sample(n_s)])
  obs_fold <- fold_map[as.character(case_ids)]
  preds    <- rep(NA_real_, nrow(X))

  for (f in seq_len(k_act)) {
    te <- which(obs_fold == f)
    tr <- which(obs_fold != f)
    if (length(te) == 0L || length(unique(y[tr])) < 2L) next
    preds[te] <- fit_ridge(X[tr, , drop = FALSE], y[tr],
                           X[te, , drop = FALSE], lam)
  }
  valid <- !is.na(preds)
  row_to_subj_auroc(case_ids[valid], preds[valid], y[valid])
}

# ---------------------------------------------------------------------------
# LOCO CV (leave-one-cluster-out; reference estimator)
# Returns list(auc, roc) for DeLong CI computation
# ---------------------------------------------------------------------------
cv_loco <- function(X, y, case_ids) {
  lam      <- fixed_lambda(length(unique(case_ids)))
  subj_ids <- unique(case_ids)
  subj_p   <- setNames(numeric(length(subj_ids)), subj_ids)

  for (id in subj_ids) {
    te <- which(case_ids == id)
    tr <- which(case_ids != id)
    if (length(unique(y[tr])) < 2L) { subj_p[as.character(id)] <- 0.5; next }
    subj_p[as.character(id)] <- mean(
      fit_ridge(X[tr, , drop = FALSE], y[tr],
                X[te, , drop = FALSE], lam))
  }

  subj_y <- tapply(y, case_ids, function(x) x[1L])
  subj_y <- subj_y[names(subj_p)]
  if (length(unique(subj_y)) < 2L) return(list(auc = NA_real_, roc = NULL))
  roc_obj <- pROC::roc(as.numeric(subj_y), as.numeric(subj_p), quiet = TRUE)
  list(auc = as.numeric(pROC::auc(roc_obj)), roc = roc_obj)
}

# =============================================================================
# 5. RUN ONE PREDICTOR CONFIGURATION
# =============================================================================

run_config <- function(dat, outcome_var, pred_cols, config_label,
                       n_reps = N_REPS, k = K_FOLDS) {

  cat(sprintf("\n  --- %s ---\n", config_label))

  X        <- as.matrix(dat[, pred_cols, drop = FALSE])
  y        <- dat[[outcome_var]]
  case_ids <- dat$Case_ID

  # Naive CV: n_reps replications
  cat("    Naive CV...")
  naive_aucs <- sapply(seq_len(n_reps), function(r) {
    cv_naive(X, y, case_ids, k = k, seed = r)
  })
  cat(" done\n")

  # Subject-level CV: n_reps replications
  cat("    Subject CV...")
  subj_aucs <- sapply(seq_len(n_reps), function(r) {
    cv_subject(X, y, case_ids, k = k, seed = r)
  })
  cat(" done\n")

  # LOCO CV: single run (deterministic)
  cat("    LOCO CV...")
  loco_res <- cv_loco(X, y, case_ids)
  loco_auc <- loco_res$auc
  loco_ci  <- if (!is.null(loco_res$roc))
                as.numeric(pROC::ci.auc(loco_res$roc, method = "delong"))
              else c(NA, NA, NA)
  cat(" done\n")

  # Optimism = CV estimate - LOCO
  naive_opt  <- naive_aucs - loco_auc
  subj_dev   <- subj_aucs  - loco_auc

  list(
    config          = config_label,
    naive_aucs      = naive_aucs,
    subj_aucs       = subj_aucs,
    loco_auc        = loco_auc,
    loco_ci         = loco_ci,
    naive_opt       = naive_opt,
    subj_dev        = subj_dev,
    mean_naive      = mean(naive_aucs, na.rm = TRUE),
    mean_subj       = mean(subj_aucs,  na.rm = TRUE),
    mean_naive_opt  = mean(naive_opt,  na.rm = TRUE),
    mean_subj_dev   = mean(subj_dev,   na.rm = TRUE),
    pi_naive        = quantile(naive_aucs, c(0.025, 0.975), na.rm = TRUE),
    pi_subj         = quantile(subj_aucs,  c(0.025, 0.975), na.rm = TRUE),
    pi_naive_opt    = quantile(naive_opt,  c(0.025, 0.975), na.rm = TRUE),
    pi_subj_dev     = quantile(subj_dev,   c(0.025, 0.975), na.rm = TRUE)
  )
}

# =============================================================================
# 6. RUN ALL THREE CONFIGURATIONS FOR BOTH OUTCOMES
# =============================================================================

# Predictor sets for each arm
fixed_cols <- fixed_vars
long_cols  <- long_vars
mixed_cols <- c(fixed_vars, long_vars)

run_all <- function(dat, outcome_var, outcome_label) {

  cat(sprintf("\n\n========== %s ==========\n", outcome_label))

  res_A <- run_config(dat, outcome_var, fixed_cols,
                      "Arm A: Fixed-Only (BWGT, GA, AP5, race/ethnicity)")
  res_B <- run_config(dat, outcome_var, long_cols,
                      "Arm B: Longitudinal-Only (7 oxygenation variables)")
  res_C <- run_config(dat, outcome_var, mixed_cols,
                      "Arm C: Mixed (fixed + longitudinal)")

  # ---------------------------------------------------------------------------
  # Print formatted results table
  # ---------------------------------------------------------------------------
  cat(sprintf("\n\n--- RESULTS TABLE: %s ---\n", outcome_label))
  cat(sprintf("%-12s  %-10s  %-8s  %-18s  %-10s  %-10s  %-18s\n",
              "Config", "Strategy", "AUROC", "95% PI/CI", "Optimism",
              "(vs LOCO)", "95% PI"))
  cat(paste(rep("-", 105), collapse = ""), "\n")

  for (res in list(res_A, res_B, res_C)) {
    lbl <- substr(res$config, 1, 8)

    # Naive CV row
    cat(sprintf("%-12s  %-10s  %8.3f  [%5.3f, %5.3f]   %+10.3f  %-10s  [%+.3f, %+.3f]\n",
                lbl, "Naive 10f",
                res$mean_naive, res$pi_naive[1], res$pi_naive[2],
                res$mean_naive_opt, "",
                res$pi_naive_opt[1], res$pi_naive_opt[2]))

    # Subject CV row
    cat(sprintf("%-12s  %-10s  %8.3f  [%5.3f, %5.3f]   %+10.3f  %-10s  [%+.3f, %+.3f]\n",
                "", "Subject 10f",
                res$mean_subj, res$pi_subj[1], res$pi_subj[2],
                res$mean_subj_dev, "",
                res$pi_subj_dev[1], res$pi_subj_dev[2]))

    # LOCO row (reference)
    cat(sprintf("%-12s  %-10s  %8.3f  [%5.3f, %5.3f]   %10s  (reference)\n",
                "", "LOCO",
                res$loco_auc, res$loco_ci[1], res$loco_ci[3], "---"))
    cat(paste(rep("-", 105), collapse = ""), "\n")
  }

  # ---------------------------------------------------------------------------
  # Additivity test: Delta_naive,mixed vs Delta_naive,A + Delta_naive,B
  # ---------------------------------------------------------------------------
  cat(sprintf("\n--- ADDITIVITY CHECK: %s ---\n", outcome_label))
  delta_A   <- res_A$mean_naive_opt
  delta_B   <- res_B$mean_naive_opt
  delta_C   <- res_C$mean_naive_opt
  additive  <- delta_A + delta_B
  departure <- delta_C - additive

  cat(sprintf("  Delta_naive (fixed-only):       %+.4f\n", delta_A))
  cat(sprintf("  Delta_naive (longitudinal-only): %+.4f\n", delta_B))
  cat(sprintf("  Sum (additive prediction):       %+.4f\n", additive))
  cat(sprintf("  Delta_naive (mixed):             %+.4f\n", delta_C))
  cat(sprintf("  Departure from additivity:       %+.4f  (%s)\n",
              departure,
              ifelse(abs(departure) < 0.01, "consistent with additivity",
                     ifelse(departure > 0, "super-additive", "sub-additive"))))

  # Also report subject-level additivity
  delta_A_s  <- res_A$mean_subj_dev
  delta_B_s  <- res_B$mean_subj_dev
  delta_C_s  <- res_C$mean_subj_dev
  add_s      <- delta_A_s + delta_B_s
  depart_s   <- delta_C_s - add_s
  cat(sprintf("\n  Subject-level deviation check:\n"))
  cat(sprintf("  Delta_subject (fixed-only):       %+.4f\n", delta_A_s))
  cat(sprintf("  Delta_subject (longitudinal-only): %+.4f\n", delta_B_s))
  cat(sprintf("  Sum (additive prediction):         %+.4f\n", add_s))
  cat(sprintf("  Delta_subject (mixed):             %+.4f\n", delta_C_s))
  cat(sprintf("  Departure from additivity:         %+.4f\n", depart_s))

  invisible(list(A = res_A, B = res_B, C = res_C))
}

# =============================================================================
# 7. EXECUTE
# =============================================================================

set.seed(20260505L)  # master seed matching simulation study

cat("\nRunning primary outcome analysis (any ROP)...\n")
res_any <- run_all(dat, "any_rop",
                   "PRIMARY: Any ROP (n=101, 48 events, 47.5%)")

cat("\nRunning secondary outcome analysis (severe ROP)...\n")
res_sev <- run_all(dat, "severe_rop",
                   "SECONDARY: Severe ROP (n=101, 15 events, 14.9%)")

# =============================================================================
# 8. SENSITIVITY: LOCO equivalence (subject-level 10f vs LOCO)
# =============================================================================

cat("\n========== LOCO vs SUBJECT 10-FOLD EQUIVALENCE (any ROP) ==========\n")
cat("(Replicates finding from companion paper; expected near-equivalence)\n\n")

for (res in list(res_any$A, res_any$B, res_any$C)) {
  diff <- res$mean_subj - res$loco_auc
  cat(sprintf("%-50s  Subject 10f - LOCO = %+.4f\n",
              res$config, diff))
}

cat(sprintf("\n\nAnalysis complete: %s\n", Sys.time()))
cat(sprintf("Output saved to: %s\n", output_file))
